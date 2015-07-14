__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob

def plotError(output, style, label):
    errMap = output['errMap']
    pyplot.plot(list(errMap[output['maxIter']].values()), style, label=label)

def plotIter(output, style, label):
    iterMap = output['iterMap']
    pyplot.plot(list(iterMap.values()), style, label=label)

def plotComparison(plotType):
    for file in glob.glob('*.dump'):
        if plotType == 'Error': pyplot.yscale('log')
        print(file)
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
        if plotType == 'Error':
            plotError(dumpObj['output'], 's-', 'err_%s' % file)
        elif plotType == 'Iter':
            plotIter(dumpObj['output'], '^-', 'iter_%s' % file)

def getSpeedUps(fileNames, sdcTime):
    speedUps = []
    for file in glob.glob(fileNames):
        print(file)
        with open(file, 'rb') as input:
            timeMeasure = pickle.load(input)['output']['timeMeasure']
            speedUp = sdcTime/timeMeasure
            print("timeMeasure: %g speedUp: %g " %(timeMeasure, speedUp))
            speedUps += [speedUp]
    return speedUps

def plotSpeedUp():
    with open('sdc.dump', 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj['input']
        sdcTime = dumpObj['output']['timeMeasure']

    speedUps = getSpeedUps('parareal_mpi_nproc*.dump', sdcTime)
    pyplot.plot(list(range(2, len(speedUps)+2)), speedUps, 's-', label='parareal_mpi')

    speedUps = getSpeedUps('parareal_hybrid_nproc*.dump', sdcTime)
    pyplot.plot(list(range(2, len(speedUps)+2)), speedUps, 's-', label='parareal_hybrid')

    pyplot.title('num_nodes: %d num_nodes_coarse: %d spatial_dofs: %d spatial_dofs_coarse: %d abs_res_tol: %g' %
                 (input.num_nodes, input.num_nodes_coarse, input.spatial_dofs, input.spatial_dofs_coarse,
                  input.abs_res_tol))

    pyplot.xlabel('Number of Processors')
    pyplot.ylabel('Speedup')

# plotSpeedUp()
plotComparison("Error")

pyplot.legend()
pyplot.show()