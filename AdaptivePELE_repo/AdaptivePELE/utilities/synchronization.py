import os
import time
import glob
import fcntl
import errno
from AdaptivePELE.utilities import utilities

try:
    ProcessLookupError
except NameError:
    ProcessLookupError = OSError


class ProcessesManager:
    """
        Object that sinchronizes multiple adaptivePELE instances, designed to
        be able to use multiple nodes of a gpu cluster
    """
    RUNNING = "RUNNING"
    WAITING = "WAITING"
    INIT = "INIT"

    def __init__(self, output_path, num_replicas):
        self.syncFolder =  os.path.join(os.path.abspath(output_path), "synchronization")
        utilities.makeFolder(self.syncFolder)
        self.lockFile = os.path.join(self.syncFolder,  "syncFile.lock")
        self.pid = os.getpid()
        self.nReplicas = num_replicas
        self.id = None
        self.lockInfo = {}
        self.status = self.INIT
        self.sleepTime = 0.5
        self.createLockFile()
        self.lock_available = True
        self.testLock()
        self.writeProcessInfo()
        self.initLockFile()
        self.syncStep = 0

    def __len__(self):
        # define the size of the ProcessesManager object as the number of
        # processes
        return len(self.lockInfo)

    def testLock(self):
        """
            Run a test to check if the file systems supports lockf and store it
            in the lock_available attribute
        """
        file_lock = open(self.lockFile, "r+")
        try:
            fcntl.lockf(file_lock, fcntl.LOCK_EX)
        except (IOError, OSError) as exc:
            if exc.errno != errno.ENOLCK:
                raise
            # if only one replica is running we don't need synchronization,
            # so we can still run even if we are using a filesystem that
            # does not support lockf
            if self.nReplicas > 1:
                raise ValueError("There was a problem allocating a lock, this usually happens when the filesystem used does not support lockf, such as NFS")
            else:
                self.lock_available = False
                return
        self.lock_available = True

    def writeProcessInfo(self):
        """
            Write the information of the running process to create the sync file
        """
        with open(os.path.join(self.syncFolder, "%d.proc" % self.pid), "w") as fw:
            fw.write("%d\n" % self.pid)

    def createLockFile(self):
        """
            Create the lock file
        """
        try:
            # whith the flags set, if the file exists the open will fail
            fd = os.open(self.lockFile, os.O_CREAT | os.O_EXCL)
            os.close(fd)
        except OSError:
            return
        file_lock = open(self.lockFile, "w")
        file_lock.write("0\n")
        file_lock.close()

    def initLockFile(self):
        """
            Initialize and write the information for the current process

        """
        while True:
            processes = glob.glob(os.path.join(self.syncFolder, "*.proc"))
            processes.sort()
            if len(processes) > self.nReplicas:
                raise utilities.ImproperParameterValueException("More processors files than replicas found, this could be due to wrong number of replicas chosen in the control file or files remaining from previous that were not clean properly")
            if len(processes) != self.nReplicas:
                time.sleep(self.sleepTime)
                continue
            # only reach this block if all processes have written their own
            # files
            for i, process_file in enumerate(processes):
                process = int(os.path.splitext(os.path.split(process_file)[1])[0])
                self.lockInfo[process] = (i, set([self.status]))
                if process == self.pid:
                    self.id = i
            break
        file_lock = open(self.lockFile, "r+")
        while True:
            if self.lock_available:
                # if the filesystem does not support locks but only one replica
                # is running we will use no locks
                fcntl.lockf(file_lock, fcntl.LOCK_EX)

            if self.isMaster():
                self.writeLockInfo(file_lock)
                if self.lock_available:
                    fcntl.lockf(file_lock, fcntl.LOCK_UN)
                file_lock.close()
                return
            else:
                lock_info = self.getLockInfo(file_lock)
                if self.lock_available:
                    fcntl.lockf(file_lock, fcntl.LOCK_UN)
                # ensure that all process are created before continuing
                if sorted(list(lock_info)) == sorted(list(self.lockInfo)):
                    file_lock.close()
                    return

    def getLockInfo(self, file_descriptor):
        """
            Return the info stored in the lock file

            :param file_descriptor: File object of the lock file
            :type file_descriptor: file

            :returns: dict -- A dictonary containing the information of the
            different process managers initialized
        """
        info = {}
        file_descriptor.seek(0)
        for line in file_descriptor:
            if line == "0\n":
                break
            pid, id_num, label = line.rstrip().split(":")
            info[int(pid)] = (int(id_num), set(label.split(";")))
        return info

    def writeLockInfo(self, file_descriptor):
        """
            Write the lock info to the lock file

            :param file_descriptor: File object of the lock file
            :type file_descriptor: file
        """
        file_descriptor.seek(0)
        for pid, data in self.lockInfo.items():
            id_num, status_set = data
            file_descriptor.write("%d:%d:%s\n" % (pid, id_num, ";".join(status_set)))
        file_descriptor.truncate()

    def isMaster(self):
        """
            Return wether the current process is the master process

            :returns: bool -- Whether the current process is master
        """
        return self.id == 0

    def setStatus(self, status):
        """
            Set the current status of the process

            :param status: Status of the process (INIT, WAITING or RUNNING)
            :type status: str
        """
        self.status = status
        file_lock = open(self.lockFile, "r+")
        if self.lock_available:
            fcntl.lockf(file_lock, fcntl.LOCK_EX)
        self.lockInfo = self.getLockInfo(file_lock)
        self.lockInfo[self.pid][1].add(self.status)
        self.writeLockInfo(file_lock)
        if self.lock_available:
            fcntl.lockf(file_lock, fcntl.LOCK_UN)
        file_lock.close()

    def getStatus(self):
        """
            Return the current status of the process

            :returns: str -- Status of the process
        """
        return self.status

    def isSynchronized(self, status):
        """
            Return wether all processes are synchronized, that is, they all
            have the same status

            :param status: Status of the process (INIT, WAITING or RUNNING)
            :type status: str

            :returns: bool -- Whether all processes are synchronized
        """
        assert len(self.lockInfo) > 0, "No processes found in lockInfo!!!"
        for pid in self.lockInfo:
            if status not in self.lockInfo[pid][1]:
                return False
        return True

    def synchronize(self, status):
        """
            Create a barrier-like situation to wait for all processes to finish
        """
        while True:
            if self.isSynchronized(status):
                return
            file_lock = open(self.lockFile, "r+")
            if self.lock_available:
                fcntl.lockf(file_lock, fcntl.LOCK_EX)
            self.lockInfo = self.getLockInfo(file_lock)
            if self.lock_available:
                fcntl.lockf(file_lock, fcntl.LOCK_UN)
            file_lock.close()

    def allRunning(self):
        """
            Check if all processes are still running
        """
        for pid in self.lockInfo:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                utilities.print_unbuffered("Process %d not found!!!" % pid)
                return False
        return True

    def writeEquilibrationStructures(self, path, structures):
        """
            Write the equilibration structures for the current replica

            :param path: Path where to write the structures
            :type path: str
            :param structures: Filename with the structures
            :type structures: list
        """
        outfile = os.path.join(path, "structures_equilibration_%d.txt" % self.id)
        with open(outfile, "w") as fw:
            fw.write(",".join(structures))

    def readEquilibrationStructures(self, path):
        """
            Read the equilibration structures for all replicas

            :param path: Path from where to read the structures
            :type path: str
        """
        files = glob.glob(os.path.join(path, "structures_equilibration_*.txt"))
        assert len(files) == self.__len__(), "Missing files for some of the replicas"
        structures = []
        for i in range(self.__len__()):
            with open(os.path.join(path, "structures_equilibration_%d.txt" % i)) as fr:
                structure_partial = fr.read().rstrip().split(",")
            if structure_partial != [""]:
                structures.extend(structure_partial)
        return structures

    def getStructureListPerReplica(self, initialStructures, trajsPerReplica):
        """
            Filter the list of initial structures to select only the ones
            corresponding to the current replica

            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: list of str
            :param trajsPerReplica: Number of trajectories that each replica has to calculate
            :type trajsPerReplica: int

            :returns: list -- List with a tuple containing the initial structures of the replica and their indices
        """
        n_structures = len(initialStructures)
        end = min(n_structures, (self.id+1)*trajsPerReplica)
        return [(i, initialStructures[i]) for i in range(self.id*trajsPerReplica, end)]

    def getBarrierName(self):
        """
            Create a unique status name so that every time we synchronize we do
            it under a different name
        """
        self.syncStep += 1
        return "%s-%d" % (self.WAITING, self.syncStep)

    def barrier(self):
        """
            Create a barrier
        """
        status = self.getBarrierName()
        self.setStatus(status)
        self.synchronize(status)
