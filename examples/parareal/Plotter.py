__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob
import re

styles = ['s-', '*-', '^-', 'o-']

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
    i = 0
    for type in ['sdc', 'parareal_mpi_nproc*', 'parareal_hybrid_nproc*', 'pfasst_nproc*']:
        for file in glob.glob('%s.dump' % type):
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
                plotError(dumpObj.output, styles[i % len(styles)], file)
            elif plotType == 'Iter':
                plotIter(dumpObj, styles[i % len(styles)], file)
            elif plotType == 'LastError':
                plotLastError(dumpObj.output, styles[i % len(styles)], file)
        i += 1


def getSpeedUps(fileNames, sdcTime):
    speedUps = []
    efficiency = []
    nprocs = []
    for file in glob.glob(fileNames):
        print(file)
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
            timeMeasure = dumpObj.output.timeMeasure
        speedUp = sdcTime/timeMeasure
        print("timeMeasure: %g speedUp: %g " %(timeMeasure, speedUp))
        speedUps += [speedUp]
        nprocs += [dumpObj.nproc]
        efficiency += [speedUp/dumpObj.nproc]
    return nprocs, speedUps, efficiency


def plotSpeedUp(efficiency):
    with open('sdc.dump', 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj.input
        sdcTime = dumpObj.output.timeMeasure

    i = 0
    for type in ['parareal_mpi', 'parareal_hybrid', 'pfasst']:
        nprocs, speedUps, efficiency  = getSpeedUps('%s_nproc*.dump' % type, sdcTime)
        if len(nprocs) > 0:
            if efficiency: plotData = efficiency
            else: plotData = speedUps
            pyplot.plot(nprocs, plotData, styles[i % len(styles)], label=type)
            i+=1

    pyplot.title('num_nodes: %d num_nodes_coarse: %d spatial_dofs: %d spatial_dofs_coarse: %d abs_res_tol: %g' %
                 (input.num_nodes, input.num_nodes_coarse, input.spatial_dofs, input.spatial_dofs_coarse,
                  input.abs_res_tol))

    pyplot.xlabel('Anzahl Prozessoren')
    if efficiency: pyplot.ylabel('Efficiency')
    else: pyplot.ylabel('Speedup')

pyplot.figure()
plotSpeedUp(True)
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