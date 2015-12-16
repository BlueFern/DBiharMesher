import argparse
import math

# The output from this script can be redirected to a file.

if __name__ == "__main__":

    argParser = argparse.ArgumentParser(
        description='Generate output containing N commands to be used by batecher with -np N where N is the number of procs.')
    argParser.add_argument('numProcs', type=int, help='Number of processors to be used with batcher.')
    argParser.add_argument('numSteps', type=int, help='Number of steps to be divided between the specified number of processors.')
    args = argParser.parse_args()

    # The floor value for the number of steps per processor.
    stepPerProc = int(math.floor(args.numSteps / float(args.numProcs)))

    for proc in range(args.numProcs):

        start = stepPerProc * proc
        adjustedStep = stepPerProc
        if proc == args.numProcs - 1:
            adjustedStep = args.numSteps - start

        end = start + adjustedStep - 1

        print 'python /hpc/scratch/bloodflo/Simulations/Flat/c4080f/CopyAttributes.py ' + str(start) + ' ' + str(end)
