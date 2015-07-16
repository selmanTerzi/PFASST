__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob
import re


def plotError(output, style, label):
    errMap = output.errMap
    pyplot.plot(list(errMap[output.maxIter].values()), style, label=label)


def plotIter(dumpObj, style, label):
    output = dumpObj.output
    iterMap = output.iterMap
    nproc = dumpObj.nproc
    iterValues = list(iterMap.values())
    color = 0
    if nproc > 1:
        maxStep = output.maxStep
        i = 0
        while i < maxStep:
            if i + nproc > maxStep:
                maxRange = maxStep - 1
            else:
                maxRange = i + nproc
            if color == 0:
                line, = pyplot.plot(list(range(i, maxRange)), iterValues[i:maxRange], style, label=label)
                color = line.get_color()
            else:
                pyplot.plot(list(range(i, maxRange)), iterValues[i:maxRange], style, color=color)
            i += nproc
    else:
        pyplot.plot(iterValues, style, label=label)


def plotLastError(output, style, label):
    errMap = output.errMap
    maxIter = output.maxIter
    maxStep = output.maxStep
    err = []
    for i in range(1, maxIter+1):
        err += [errMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter+1)), err, style, label=label)


def plotComparison(plotType):
    for file in glob.glob('*.dump'):
        if plotType in ['Error', 'LastError']:
            pyplot.yscale('log')
            pyplot.ylabel('Fehler')
            if plotType == 'LastError':
                pyplot.xlabel('Iterationen')
        if plotType in ['Error', 'Iter']:
            pyplot.xlabel('Zeitschritt')
            if plotType == 'Iter':
                pyplot.ylabel('Iterationen')
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
        file = file.split('.')[0]
        if plotType == 'Error':
            plotError(dumpObj.output, 's-', file)
        elif plotType == 'Iter':
            plotIter(dumpObj, '^-', file)
        elif plotType == 'LastError':
            plotLastError(dumpObj.output, '^-', file)


def getSpeedUps(fileNames, sdcTime):
    speedUps = []
    nprocs = []
    for file in glob.glob(fileNames):
        print(file)
        nproc = int(re.findall("nproc(\d*)", file)[0])
        with open(file, 'rb') as input:
            timeMeasure = pickle.load(input).output.timeMeasure
        speedUp = sdcTime/timeMeasure
        print("timeMeasure: %g speedUp: %g " %(timeMeasure, speedUp))
        speedUps += [speedUp]
        nprocs += [nproc]
    return speedUps, nprocs


def plotSpeedUp():
    with open('sdc.dump', 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj.input
        sdcTime = dumpObj.output.timeMeasure

    speedUps, nprocs  = getSpeedUps('parareal_mpi_nproc*.dump', sdcTime)
    if len(speedUps) > 0:
        pyplot.plot(nprocs, speedUps, 's-', label='parareal_mpi')

    speedUps, nprocs = getSpeedUps('parareal_hybrid_nproc*.dump', sdcTime)
    if len(speedUps) > 0:
        pyplot.plot(nprocs, speedUps, 's-', label='parareal_hybrid')

    speedUps, nprocs = getSpeedUps('pfasst_nproc*.dump', sdcTime)
    if len(speedUps) > 0:
        pyplot.plot(nprocs, speedUps, 's-', label='pfasst')

    pyplot.title('num_nodes: %d num_nodes_coarse: %d spatial_dofs: %d spatial_dofs_coarse: %d abs_res_tol: %g' %
                 (input.num_nodes, input.num_nodes_coarse, input.spatial_dofs, input.spatial_dofs_coarse,
                  input.abs_res_tol))

    pyplot.xlabel('Anzahl Prozessoren')
    pyplot.ylabel('Speedup')

pyplot.figure()
plotSpeedUp()
pyplot.legend()

pyplot.figure()
plotComparison("LastError")
pyplot.legend()

pyplot.figure()
plotComparison("Error")
pyplot.legend()

pyplot.figure()
plotComparison("Iter")
pyplot.legend()

pyplot.show()