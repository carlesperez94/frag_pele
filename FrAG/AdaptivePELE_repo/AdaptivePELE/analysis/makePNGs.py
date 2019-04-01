"""
    Create two PNG images with RMSD-MC steps and BE-RMSD for a set of different
    adaptive runs
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import subprocess
import glob
import argparse
from AdaptivePELE.analysis import plotAdaptive


def parseArguments():
    """
        Parse command line arguments
    """
    desc = "It makes two PNG images with RMSD-steps and BE-RMSD for a set of different adaptive runs\n\n"\
           "Instructions:\n"\
           "-------------\n"\
           "Keys in \"folders\", \"titles\" and \"outputFilenames\" dictionaries should match (could be avoided with simple refactor...)\n"\
           "Keys are: \"#steps_#processors\", #steps is important for the plots"\
           "Change params to match simulation.\n"
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    folders = {"4_512": "PRprog_4_512"}  # , "4_64":"3ptb_4_64", "4_128":"3ptb_4_128"}

    subfoldersWildcard = "inversely_*"
    subfoldersWildcard = "be_epsilon_*"

    titles = {"4_512": "n=32, 4 steps, %s"}  # , "4_64":"n=64, 4 steps, %s", "4_128":"n=128, 4 steps, %s"}

    outputFilenames = {"4_512": "512_4_%s"}  # , "4_64":"64_4_%s", "4_128":"128_4_%s"}

    parameters = {"stepsCol": 2,
                  "RMSDCol": 5,
                  "BECol": 6,
                  "reportFilename": "report_"}

    gplFolder = "/gpfs/scratch/bsc72/bsc72755/adaptiveSampling/simulation"
    tmpFolder = "/tmp"

    tmpPlotFile = os.path.join(tmpFolder, "tmp.gpl")

    gnuplot = "$SCRATCH/software/gnuplot/bin/gnuplot"

    def buildGnuplotString(title, outputFilename, params):
        gnuplotFileStringContent = """\
        set term png\n\
        set title "%(plotTitle)s"\n\
        set output "rmsd_steps_%(outputFilename)s.png"\n\
        %(rmsdStepsPringString)s\n\
        \n\
        set output "be_rmsd_%(outputFilename)s.png\n\
        %(beRmsdPrintString)s\n
        """

        stepsPerRun = params["stepsPerRun"]
        stepsCol = params["stepsCol"]
        RMSDCol = params["RMSDCol"]
        BECol = params["BECol"]
        reportFilename = params["reportFilename"]

        rmsdStepsPrintString = plotAdaptive.generatePrintString(stepsPerRun, stepsCol, RMSDCol, reportFilename, "PRINT_RMSD_STEPS")
        beRmsdPrintString = plotAdaptive.generatePrintString(stepsPerRun, RMSDCol, BECol, reportFilename, "PRINT_BE_RMSD")
        dictionary = {"plotTitle": title,
                      "outputFilename": outputFilename,
                      "rmsdStepsPringString": rmsdStepsPrintString,
                      "beRmsdPrintString": beRmsdPrintString}

        return gnuplotFileStringContent % dictionary

    parseArguments()

    for key, folder in folders.items():
        print("Folder: ", folder)
        try:
            os.chdir(folder)
        except OSError:
            continue

        parameters["stepsPerRun"] = int(key.split("_")[0])

        subfolders = glob.glob(subfoldersWildcard)
        print(subfolders)
        for subfolder in subfolders:
            os.chdir(subfolder)

            gnuplotFileContent = buildGnuplotString(titles[key] % subfolder, outputFilenames[key] % subfolder, parameters)

            with open(tmpPlotFile, "w") as f:
                f.write(gnuplotFileContent)

            proc = subprocess.Popen("%s %s" % (gnuplot, tmpPlotFile), stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            if out:
                print(out)
            if err:
                print(err)
            os.chdir("..")

        os.chdir("..")
