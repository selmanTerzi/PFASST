__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob
from ProcessStarter import RunTypes
from numpy import log10

styles = ['s-', '*-', '^-', 'o-', '+-']
stylecounter = 0


def nextStyle():
    global stylecounter
    s = styles[stylecounter % len(styles)]
    stylecounter += 1
    return s


class PlotTypes:
    AllIterationResiduals = 1
    AllIterationErrors = 2
    Iter = 3
    LastStepError = 4
    LastStepResidual = 5
    LastIterationError = 6
    SpeedUp = 7
    Efficiency = 8
    ConvergenceData = 9


def plotLastIterationError(dumpObj, style, label):
    errMap = dumpObj.result.errMap
    pyplot.plot(list(errMap[dumpObj.result.maxIter].values()), style, label=label)


def plotAllIterations(data, style, label):
    firstPlot = True
    for i in data.keys():
        if firstPlot:
            line, = pyplot.plot(list(data[i].values()), style, label=label)
            color = line.get_color()
            firstPlot = False
        else:
            pyplot.plot(list(data[i].values()), style, color=color)


def plotIterations(dumpObj, style, label):
    result = dumpObj.result
    iterMap = result.iterMap
    nproc = dumpObj.nproc
    iterValues = list(iterMap.values())
    color = 0
    if nproc > 1:
        maxStep = result.maxStep
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
    errMap = dumpObj.result.errMap
    maxIter = dumpObj.result.maxIter
    maxStep = dumpObj.result.maxStep
    err = []
    for i in range(1, maxIter+1):
        err += [errMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter+1)), err, style, label=label)

def plotLastStepResidual(dumpObj, style, label):
    resMap = dumpObj.result.resMap
    maxIter = dumpObj.result.maxIter
    maxStep = dumpObj.result.maxStep
    err = []
    for i in range(1, maxIter+1):
        err += [resMap[i][maxStep]]
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


def plotComparison(plotType, plotData, nprocs = [0]):
    sdcTypes = [RunTypes.SDC_Coarse, RunTypes.SDC_Fine]
    for nproc in nprocs:
        for fname in processPlotDataName(set(plotData).difference(sdcTypes), nproc):
            for f in glob.glob(fname):
                plotFile(plotType, f)
    for pdata in set(plotData).intersection(sdcTypes):
        plotFile(plotType, pdata + '.dump')


plotDict = {PlotTypes.LastIterationError: plotLastIterationError,
            PlotTypes.Iter: plotIterations,
            PlotTypes.LastStepError: plotLastStepError,
            PlotTypes.LastStepResidual: plotLastStepResidual,
            PlotTypes.AllIterationErrors: plotAllIterations,
            PlotTypes.AllIterationResiduals: plotAllIterations}


def plotFile(plotType, fname):
    print(fname)
    if plotType in [PlotTypes.LastStepError,
                    PlotTypes.LastStepResidual,
                    PlotTypes.LastIterationError,
                    PlotTypes.AllIterationErrors,
                    PlotTypes.AllIterationResiduals,
                    PlotTypes.Iter]:
        if plotType != PlotTypes.Iter:
            pyplot.yscale('log')
        if plotType in [PlotTypes.AllIterationResiduals,
                        PlotTypes.LastStepResidual]:
            pyplot.ylabel('Residuum')
        elif plotType == PlotTypes.Iter:
            pyplot.ylabel('Iterationen')
        else:
            pyplot.ylabel('Fehler')
    if plotType in [PlotTypes.LastStepError,
                    PlotTypes.LastStepResidual]:
        pyplot.xlabel('Iterationen')
    else:
        pyplot.xlabel('Zeitschritt')

    with open(fname, 'rb') as input:
        dumpObj = pickle.load(input)

    fname = fname.split('.')[0]
    data = dumpObj
    if plotType == PlotTypes.AllIterationErrors:
        data = data.result.errMap
    elif plotType == PlotTypes.AllIterationResiduals:
        data = data.result.resMap
    plotDict[plotType](data, nextStyle(), fname)


def getSpeedUps(fileNames, sdcTime):
    data = {}
    for file in glob.glob(fileNames):
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
            timeMeasure = dumpObj.result.timeMeasure
        speedUp = sdcTime/timeMeasure
        efficiency = speedUp/dumpObj.nproc
        # print("timeMeasure: %g speedUp: %g " %(timeMeasure, speedUp))
        data[dumpObj.nproc] = [speedUp, efficiency]
        sortedData = sorted(data.items())
        nprocs = []
        speedUps = []
        efficiency = []
        for t in sortedData:
            nprocs += [t[0]]
            speedUps += [t[1][0]]
            efficiency += [t[1][1]]
    return nprocs, speedUps, efficiency


def plotSpeedUp(plotEfficiency, plotData):
    with open('%s.dump' % RunTypes.SDC_Fine, 'rb') as input:
        dumpObj = pickle.load(input)
        input = dumpObj.result.input
        sdcTime = dumpObj.result.timeMeasure

    for type in list(set([RunTypes.PARA_CLASSIC, RunTypes.PARA_HYBRID_FULL,
                          RunTypes.PARA_HYBRID_PARTIAL, RunTypes.PFASST]) & set(plotData)):
        nprocs, speedUps, efficiency = getSpeedUps('%s_nproc*.dump' % type, sdcTime)
        if len(nprocs) > 0:
            if plotEfficiency:
                plotData = efficiency
            else:
                plotData = speedUps
            pyplot.plot(nprocs, plotData, nextStyle(), label=type)

    pyplot.title('num_nodes: %d num_nodes_coarse: %d spatial_dofs: %d spatial_dofs_coarse: %d abs_res_tol: %g' %
                 (input.num_nodes, input.num_nodes_coarse, input.spatial_dofs, input.spatial_dofs_coarse,
                  input.abs_res_tol))

    pyplot.xlabel('Anzahl Prozessoren')
    if plotEfficiency:
        pyplot.ylabel('Effizienz')
    else:
        pyplot.ylabel('Speedup')


def plotConvergenceData(plotData):
    for file in glob.glob('orderPlot_numNodes[1-9].pkl'):
        with open(file, 'rb') as f:
            [errDict, dtArr, input] = pickle.load(f)

        for pdata in plotData:
            if pdata not in errDict.keys():
                continue
            errList = errDict[pdata]
            order = abs((log10(errList[1])-log10(errList[0]))/(log10(dtArr[1])-log10(dtArr[0])))

            pyplot.plot(dtArr, errList, nextStyle(), label='%s_num_nodes_%d' % (pdata, input.num_nodes))

    pyplot.xscale('log')
    pyplot.yscale('log')
    ax = pyplot.gca()
    ax.invert_xaxis()

    pyplot.xlabel('dt')
    pyplot.ylabel('Fehler')


def plot(plotTypes, plotData, nprocs=[0]):
    for pType in plotTypes:
        pyplot.figure()
        if pType in [PlotTypes.SpeedUp, PlotTypes.Efficiency]:
            plotSpeedUp(pType == PlotTypes.Efficiency, plotData)
        elif pType == PlotTypes.ConvergenceData:
            plotConvergenceData(plotData)
        else:
            plotComparison(pType, plotData, nprocs)
        pyplot.legend()
    pyplot.show()
