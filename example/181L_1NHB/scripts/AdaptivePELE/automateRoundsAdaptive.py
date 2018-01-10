from simulation import simulationrunner
import adaptiveSampling
import argparse
import os


def automateSimulation(args):
    controlFile = args.controlFile
    numSimulations = args.numSimulations
    nProcessors = args.nProcessors
    nSteps = args.nSteps
    simulationName = args.simulationName
    simulationParameters = simulationrunner.SimulationParameters()
    simulationParameters.templetizedControlFile = controlFile
    simulationRunner = simulationrunner.SimulationRunner(simulationParameters)
    epochs = args.epochs
    if epochs:
        rangeOfEpochs = epochs
    else:
        rangeOfEpochs = range(1, numSimulations+1)
    print "rangeOfEpochs", rangeOfEpochs

    for i in rangeOfEpochs:
        controlFileDictionary = {"SEED": "%d%d%d", "OUTPUTPATH": "%s_%d"}
        SEED_i = int(controlFileDictionary["SEED"] % (i, nProcessors, nSteps))
        controlFileDictionary["SEED"] = SEED_i
        outputPath_i = controlFileDictionary["OUTPUTPATH"] % (simulationName, i)
        controlFileDictionary["OUTPUTPATH"] = outputPath_i
        controlFileName = "tmp_%s_controlfile_%s_%d.conf" % (os.path.splitext(controlFile)[0], simulationName, i)
        controlFileName = controlFileName.replace("/","_")
        simulationRunner.makeWorkingControlFile(controlFileName, controlFileDictionary)
        print "Starting simulation %d" % i
        adaptiveSampling.main(controlFileName)


def parseArguments():
    parser = argparse.ArgumentParser(description="Automate the process "
                                     "of repeating simulations")
    parser.add_argument('controlFile', type=str)
    parser.add_argument('numSimulations', type=int)
    parser.add_argument('nProcessors', type=int)
    parser.add_argument('nSteps', type=int)
    parser.add_argument('simulationName', type=str)
    parser.add_argument("-e", "--epochs", nargs='*', type=int, help="Epochs to run")
    args = parser.parse_args()
    return args


def main():
    args = parseArguments()
    automateSimulation(args)

if __name__ == "__main__":
    main()
