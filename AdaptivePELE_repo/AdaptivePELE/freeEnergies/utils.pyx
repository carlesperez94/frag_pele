import numpy as np
from io import open
cimport cython
cimport numpy as np
from libc.math cimport sqrt
from libc.stdlib cimport rand, RAND_MAX


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] sum_rows(double[:,:] C):
    cdef Py_ssize_t sx, sy, i, j
    cdef double tmp_sum
    sx = C.shape[0]
    sy = C.shape[1]
    cdef double[:] x = np.empty((sy))
    for i in range(sx):
        tmp_sum = 0.0
        for j in range(sx):
            tmp_sum += C[i, j]
        x[i] = tmp_sum
    return x


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def buildRevTransitionMatrix(double[:,:] C):
    """
        #Implemented as Prinz paper
    """
    cdef double[:,:] X, T
    cdef Py_ssize_t sx, sy, k, ky
    sx = C.shape[0]
    sy = C.shape[1]
    X = np.empty((sx, sy))
    T = np.empty((sx, sy))
    for k in range(sx):
        for ky in range(sy):
            X[k, ky]=C[k,ky]+C[ky,k]
            T[k, ky]= 0
    cdef double[:] x, c
    x = sum_rows(X)
    c = sum_rows(C)
    cdef int iterations = 1000
    cdef int it, i, j
    cdef double a, b, z, tmp, tmp2, tmp3

    for it in range(iterations):
        # if it != 0 and it % 10 == 0: print it
        for i in range(sx):
            X[i,i] = C[i,i] * (x[i] - X[i,i]) / (c[i] - C[i,i])
        x = sum_rows(X)
        for i in range(sx - 1):
            for j in range(i + 1, sx):
                a = c[i] - C[i,j] + c[j] - C[j,i]
                b = c[i]*(x[j]-X[i,j]) + c[j]*(x[i] - X[i,j]) -(C[i,j] + C[j,i])*(x[i] + x[j] - 2*X[i,j])
                z = -(C[i,j] + C[j,i])*(x[i] - X[i,j])*(x[j] - X[i,j])
                tmp = b**2 - 4*a*z
                tmp2 = 2*a
                tmp3 = sqrt(tmp)
                X[j,i] = (-b+tmp3)/tmp2
                X[i,j] = X[j,i]
        x = sum_rows(X)

    for i in range(sx):
        for j in range(sx):
            T[i,j] = X[i,j]/x[i]
    return np.array(T)


@cython.boundscheck(False)
@cython.wraparound(False)
def isAlphaCarbon(basestring string, bint writeCA):
    cdef basestring CA = u"CA"
    cdef basestring C = u"C"

    return writeCA and string[12:16].strip() == CA and string[76:80].strip() == C


@cython.boundscheck(False)
@cython.wraparound(False)
def extraAtomCheck(basestring line, dict extraAtoms):
    cdef bint result
    cdef basestring resname = line[17:20].strip()
    cdef basestring atomname = line[12:16].strip()
    if resname not in extraAtoms:
        return False
    cdef basestring extra_atom = extraAtoms[resname]
    result = atomname == extra_atom
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def isSidechain(basestring string, bint writeSide, list sidechains):
    cdef basestring index = (string[6:11].strip())
    if not writeSide or not index:
        return False
    else:
        return int(index) in sidechains

@cython.boundscheck(False)
@cython.wraparound(False)
def is_model(basestring line):
    cdef basestring check = "ENDMDL"
    return check in line[:7]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_end(basestring line):
    cdef basestring check = "END"
    return check == line[:3]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_remark(basestring line):
    cdef basestring check = "REMARK"
    return check == line[:6]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_cryst(basestring line):
    cdef basestring check = "CRYST"
    return check == line[:5]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef calculateAutoCorrelation(list lagtimes, list dtrajs, int nclusters, int nLags):
    cdef double[:, :] autoCorr = np.zeros((nclusters, nLags))
    cdef int N = 0
    cdef int v1, v2
    cdef double[:] M = np.zeros(nLags)
    cdef long[:] traj
    cdef int lagtime
    cdef Py_ssize_t il, i, j, Nt
    cdef double[:] C = np.zeros(nclusters)
    cdef double[:] mean = np.zeros(nclusters)
    cdef double[:] var = np.zeros(nclusters)
    cdef double N_f, var_tmp
    cdef int maxLag = max(lagtimes)
    for traj in dtrajs:
        Nt = traj.size
        if Nt < maxLag:
            raise ValueError("Lagtime specified are too big for the trajectories!")
        N += Nt
        for i in range(Nt):
            C[traj[i]] += 1
    N_f = <double>(N)
    var_tmp = N_f*(N_f-1)
    for i in range(nclusters):
        mean[i] = C[i]/N_f
        var[i] = (N*C[i]-(C[i]**2))/var_tmp

    for traj in dtrajs:
        Nt = traj.size
        for il in range(nLags):
            lagtime = lagtimes[il]
            M[il] += Nt-lagtime
            for i in range(Nt-lagtime):
                for c in range(nclusters):
                    v1 = 0
                    v2 = 0
                    if c == traj[i]:
                        v1 = 1
                    if c == traj[i+lagtime]:
                        v2 = 1
                    autoCorr[c, il] += (v1-mean[c])*(v2-mean[c])

    for i in range(nclusters):
        for j in range(nLags):
            autoCorr[i, j] /= M[j]
            autoCorr[i, j] /= var[i]
    return np.array(autoCorr)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double[:] runSimulation(double[:, :] P, int steps, int startingPosition, long[:] states):
    cdef Py_ssize_t n = P.shape[0]
    cdef int position = startingPosition
    cdef Py_ssize_t step

    cdef double[:] traj = np.zeros(steps)
    traj[0] = position
    for step in range(1, steps):
        position = sample(states, P[position])
        traj[step] = position

    return traj

cpdef int sample(long[:] states, double[:] prob):
    """ Method obtained from http://keithschwarz.com/darts-dice-coins/ named roulette wheel selection"""
    cdef double rnd = rand() / (RAND_MAX + 1.0)
    return binary_search(rnd, prob)

cdef int binary_search(double rnd, double[:] prob):
    cdef Py_ssize_t l = 0
    cdef Py_ssize_t r = prob.shape[0]-1
    cdef int mid
    while l < r:
        mid =  (l+r) / 2
        if prob[mid] < rnd:
            l = mid + 1
        else:
            r = mid - 1
    if prob[l] < rnd:
        # return the first element that is larger than the random number
        return l + 1
    else:
        return l

@cython.boundscheck(False)
@cython.wraparound(False)
def contactMap(double[:,:] coords, Py_ssize_t n):
    cdef Py_ssize_t n_rows, size , row_ind, i, j, el
    cdef double tmp = 0.5*n*(n-1)
    size = <Py_ssize_t>tmp
    n_rows = coords.shape[0]
    cdef double tmp_sum, dist_x, dist_y, dist_z
    cdef double[:, :] cm = np.zeros((n_rows, size))
    for row_ind in range(n_rows):
        # row_ind corresponds to a frame
        el = 0
        for i in range(n):
            for j in range(i+1, n):
                dist_x = coords[row_ind][i*3]-coords[row_ind][j*3]
                dist_y = coords[row_ind][i*3+1]-coords[row_ind][j*3+1]
                dist_z = coords[row_ind][i*3+2]-coords[row_ind][j*3+2]
                cm[row_ind][el] = sqrt(dist_x**2+dist_y**2+dist_z**2)
                el += 1
    return np.asarray(cm)
