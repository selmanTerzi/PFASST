__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob
from ProcessStarter import RunTypes

styles = ['s-', '*-', '^-', 'o-', '+-']
stylecounter = 0

class PlotTypes:
    AllIterationResiduals = 1
    AllIterationErrors = 2
    Iter = 3
    LastStepError = 4
    LastIterationError = 5
    SpeedUp = 6
    Efficiency = 7


def plotLastIterationError(dumpObj, style, label):
    errMap = dumpObj.output.errMap
    pyplot.plot(list(errMap[dumpObj.output.maxIter].values()), style, label=label)


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


def plotLastStepError(dumpObj, style, label):
    errMap = dumpObj.output.errMap
    maxIter = dumpObj.output.maxIter
    maxStep = dumpObj.output.maxStep
    err = []
    for i in range(1, maxIter+1):
        err += [errMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter+1)), err, style, label=label)


def processPlotDataName(plotData, nproc):
    plotDataNames = []
    if nproc == 0:
        nproc = '*'
    else:
        nproc = '%d' % nproc
    for pdata in plotData:
        plotDataNames += ['%s_nproc%s.dump' % (pdata, nproc)]
    return plotDataNames

plotDict = {PlotTypes.LastIterationError: plotLastIterationError,
            PlotTypes.Iter: plotIters,
            PlotTypes.LastStepError: plotLastStepError,
            PlotTypes.AllIterationErrors: plotAllIterations,
            PlotTypes.AllIterationResiduals: plotAllIterations}


def plotComparison(plotType, plotData, nprocs = [0]):
    sdcTypes = [RunTypes.SDC_Coarse, RunTypes.SDC_Fine]
    for nproc in nprocs:
        for fname in processPlotDataName(set(plotData).difference(sdcTypes), nproc):
            for f in glob.glob(fname):
                plotFile(plotType, f)
    for pdata in set(plotData).intersection(sdcTypes):
        plotFile(plotType, pdata + '.dump')


def plotFile(plotType, fname):
    global stylecounter
    print(fname)
    if plotType in [PlotTypes.LastStepError,
                    PlotTypes.LastIterationError,
                    PlotTypes.AllIterationErrors,
                    PlotTypes.AllIterationResiduals,
                    PlotTypes.Iter]:
        if plotType != PlotTypes.Iter:
            pyplot.yscale('log')
        if plotType == PlotTypes.AllIterationResiduals:
            pyplot.ylabel('Residuum')
        elif plotType == PlotTypes.Iter:
            pyplot.ylabel('Iterationen')
        else:
            pyplot.ylabel('Fehler')
    if plotType == PlotTypes.LastStepError:
        pyplot.xlabel('Iterationen')
    else:
        pyplot.xlabel('Zeitschritt')

    with open(fname, 'rb') as input:
        dumpObj = pickle.load(input)

    fname = fname.split('.')[0]
    data = dumpObj
    if plotType == PlotTypes.AllIterationErrors:
        data = data.output.errMap
    elif plotType == PlotTypes.AllIterationResiduals:
        data = data.output.resMap
    plotDict[plotType](data, styles[stylecounter % len(styles)], fname)
    stylecounter +=1

def getSpeedUps(fileNames, sdcTime):
    data = {}
    for file in glob.glob(fileNames):
        print(file)
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
            timeMeasure = dumpObj.output.timeMeasure
        speedUp = sdcTime/timeMeasure
        efficiency = speedUp/dumpObj.nproc
        print("timeMeasure: %g speedUp: %g " %(timeMeasure, speedUp))
        data[dumpObj.nproc] = [speedUp, efficiency]
        sortedData = sorted(data.items())
        nprocs = []
        speedUps = []
        efficiency = []
        for t in sortedData:
            print(t)
            nprocs += [t[0]]
            speedUps += [t[1][0]]
            efficiency += [t[1][1]]
    return nprocs, speedUps, efficiency


def plotSpeedUp(plotEfficiency, plotData):
    with open('%s.dump' % RunTypes.SDC_Fine, 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj.input
        sdcTime = dumpObj.output.timeMeasure

    i = 0
    for type in list(set([RunTypes.PARA_CLASSIC, RunTypes.PARA_HYBRID_FULL,
                          RunTypes.PARA_HYBRID_PARTIAL, RunTypes.PFASST]) & set(plotData)):
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


def plot(plotTypes, plotData, nprocs=[0]):
    for pType in plotTypes:
        pyplot.figure()
        if pType in [PlotTypes.SpeedUp, PlotTypes.Efficiency]:
            plotSpeedUp(pType == PlotTypes.Efficiency, plotData)
        else:
            plotComparison(pType, plotData, nprocs)
        pyplot.legend()
    pyplot.savefig('out')
    pyplot.show()