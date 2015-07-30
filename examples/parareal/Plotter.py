__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob
from ConvenientProcessStarter import RunTypes

styles = ['s-', '*-', '^-', 'o-', '+-']


class PlotTypes:
    AllIterationResiduals = 1
    AllIterationErrors = 2
    Iter = 3
    LastStepError = 4
    LastIterationError = 5
    SpeedUp = 6
    Efficiency = 7


def plotLastIterationError(output, style, label):
    errMap = output.errMap
    pyplot.plot(list(errMap[output.maxIter].values()), style, label=label)


def plotAllIterations(data, style, label):
    firstPlot = True
    for i in data.keys():
        if firstPlot:
            line, = pyplot.plot(list(data[i].values()), style, label=label)
            color = line.get_color()
            firstPlot = False
        else:
            pyplot.plot(list(data[i].values()), style, color=color)


def plotIters(dumpObj, style, label):
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


def plotLastStepError(output, style, label):
    errMap = output.errMap
    maxIter = output.maxIter
    maxStep = output.maxStep
    err = []
    for i in range(1, maxIter+1):
        err += [errMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter+1)), err, style, label=label)


def getRegex(data):
    np = '*'
    if data != 0:
        if isinstance(data, list) > 0:
            np = "["
            for i in data:
                np += "%d" % i
            np += "]"
        else:
            np = "%d" % data
    return np


def processPlotDataNames(plotData, nprocs):
    plotDataNames = []
    np = getRegex(nprocs)
    for pdata in plotData:
        if pdata in [RunTypes.SDC_Fine, RunTypes.SDC_Coarse]:
            plotDataNames += ['%s.dump' % pdata]
        else:
            plotDataNames += ['%s_nproc%s.dump' % (pdata,np)]
    return plotDataNames


def plotComparison(plotType, plotData, nprocs = 0):
    i = 0
    for fname in processPlotDataNames(plotData, nprocs):
        for file in glob.glob(fname):
            if plotType in [PlotTypes.LastStepError,
                            PlotTypes.LastIterationError,
                            PlotTypes.AllIterationErrors,
                            PlotTypes.AllIterationResiduals]:
                pyplot.yscale('log')
                if plotType == PlotTypes.AllIterationResiduals:
                    pyplot.ylabel('Residuum')
                else:
                    pyplot.ylabel('Fehler')
                if plotType == PlotTypes.LastIterationError:
                    pyplot.xlabel('Iterationen')
                else:
                    pyplot.xlabel('Zeitschritt')
            elif plotType == PlotTypes.Iter:
                pyplot.xlabel('Zeitschritt')
                pyplot.ylabel('Iterationen')

            with open(file, 'rb') as input:
                dumpObj = pickle.load(input)
            file = file.split('.')[0]
            if plotType == PlotTypes.LastIterationError:
                plotLastIterationError(dumpObj.output, styles[i % len(styles)], file)
            elif plotType == PlotTypes.Iter:
                plotIters(dumpObj, styles[i % len(styles)], file)
            elif plotType == PlotTypes.LastStepError:
                plotLastStepError(dumpObj.output, styles[i % len(styles)], file)
            elif plotType == PlotTypes.AllIterationErrors:
                plotAllIterations(dumpObj.output.errMap, styles[i % len(styles)], file)
            elif plotType == PlotTypes.AllIterationResiduals:
                plotAllIterations(dumpObj.output.resMap, styles[i % len(styles)], file)
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


def plotSpeedUp(plotEfficiency, plotData):
    print(plotEfficiency)
    with open('%s.dump' % RunTypes.SDC_Fine, 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj.input
        sdcTime = dumpObj.output.timeMeasure

    i = 0
    for type in list(set([RunTypes.PARA_CLASSIC, RunTypes.PARA_HYBRID_FULL, RunTypes.PFASST]) & set(plotData)):
        nprocs, speedUps, efficiency = getSpeedUps('%s_nproc*.dump' % type, sdcTime)
        if len(nprocs) > 0:
            if plotEfficiency:
                plotData = efficiency
            else:
                plotData = speedUps
            pyplot.plot(nprocs, plotData, styles[i % len(styles)], label=type)
            i+=1

    pyplot.title('num_nodes: %d num_nodes_coarse: %d spatial_dofs: %d spatial_dofs_coarse: %d abs_res_tol: %g' %
                 (input.num_nodes, input.num_nodes_coarse, input.spatial_dofs, input.spatial_dofs_coarse,
                  input.abs_res_tol))

    pyplot.xlabel('Anzahl Prozessoren')
    if plotEfficiency:
        pyplot.ylabel('Effizienz')
    else:
        pyplot.ylabel('Speedup')


def plot(plotTypes, plotData, nprocs=0):
    for pType in plotTypes:
        pyplot.figure()
        if pType in [PlotTypes.SpeedUp, PlotTypes.Efficiency]:
            plotSpeedUp(pType == PlotTypes.Efficiency, plotData)
        else:
            plotComparison(pType, plotData, nprocs)
        pyplot.legend()
    pyplot.show()