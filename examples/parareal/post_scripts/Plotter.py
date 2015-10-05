import pickle
import glob
import numpy
import math
import re
from matplotlib import pyplot
from ProcessStarter import RunTypes
from numpy import log10


# class containing an id for each plot type provided by this plotter script
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
    SpeedUpAlpha = 10
    Timings = 11


labelDict = {RunTypes.SDC_Fine: 'feines SDC',
             RunTypes.SDC_Coarse: 'grobes SDC',
             RunTypes.PARA_CLASSIC_SERIAL: 'serielles klassisches Parareal',
             RunTypes.PARA_CLASSIC: 'klassisches Parareal',
             RunTypes.PARA_HYBRID_FULL: 'voll-hybrides Parareal',
             RunTypes.PARA_HYBRID_PARTIAL: 'teil-hybrides Parareal',
             RunTypes.PFASST: 'PFASST'}

styleDict = {RunTypes.SDC_Fine: 'o',
             RunTypes.SDC_Coarse: 'v',
             RunTypes.PARA_CLASSIC_SERIAL: 'p',
             RunTypes.PARA_CLASSIC: 's',
             RunTypes.PARA_HYBRID_FULL: '^',
             RunTypes.PARA_HYBRID_PARTIAL: '*',
             RunTypes.PFASST: '+'}

# set containing the runtypes for which the theoretical speedups can be plottet
speedUpAlphaSet = set([RunTypes.PARA_CLASSIC, RunTypes.PARA_HYBRID_FULL, RunTypes.PARA_HYBRID_PARTIAL])
# set containing the runtypes for which speedups can be plottet
speedUpSet = speedUpAlphaSet | set(RunTypes.PFASST)


def getStyle(runtype, dashed):
    if dashed:
        line = '--'
    else:
        line = '-'
    return styleDict[runtype] + line


def plotLastIterationError(dumpObj, style, label):
    errMap = dumpObj.result.errMap
    pyplot.plot(list(errMap[dumpObj.result.maxIter].values()), style, label=label)
    return 1


def plotAllIterations(data, style, label):
    firstPlot = True
    for i in data.keys():
        if firstPlot:
            line, = pyplot.plot(list(data[i].values()), style, label=label)
            color = line.get_color()
            firstPlot = False
        else:
            pyplot.plot(list(data[i].values()), style, color=color)
    return 1


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
    return 1


def plotLastStepError(dumpObj, style, label):
    errMap = dumpObj.result.errMap
    maxIter = dumpObj.result.maxIter
    maxStep = dumpObj.result.maxStep
    err = []
    for i in range(1, maxIter + 1):
        err += [errMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter + 1)), err, style, label=label)
    return 1


def plotTimings(dumpObj, style, label):
    timing = dumpObj.result.timing
    nproc = dumpObj.nproc
    x = list(range(1, nproc+1))
    pyplot.plot(x, timing.totalTimes, style, label=labelDict[dumpObj.runtype] + ' GesamtZeit')
    pyplot.plot(x, timing.interpolTimes, style, label=labelDict[dumpObj.runtype] + ' Interpolation')
    pyplot.plot(x, timing.communicationTimes, style, label=labelDict[dumpObj.runtype] + ' Kommunikation')
    return 3


def plotLastStepResidual(dumpObj, style, label):
    resMap = dumpObj.result.resMap
    maxIter = dumpObj.result.maxIter
    maxStep = dumpObj.result.maxStep
    err = []
    for i in range(1, maxIter + 1):
        err += [resMap[i][maxStep]]
    pyplot.plot(list(range(1, maxIter + 1)), err, style, label=label)
    return 1


def processPlotDataName(plotData, nproc, root='.'):
    plotDataNames = []
    if nproc == 0:
        nproc = '*'
    else:
        nproc = '%d' % nproc
    for pdata in plotData:
        plotDataNames += ['%s/%s_nproc%s.dump' % (root, pdata, nproc)]
    return plotDataNames


def plotComparison(plotType, plotData, nprocs=[0], root='.'):
    numPlots = 0
    for nproc in nprocs:
        for fname in processPlotDataName(set(plotData).difference([RunTypes.SDC_Coarse, RunTypes.SDC_Fine]), nproc, root):
            for f in glob.glob(fname):
                numPlots += plotFile(plotType, f)
    for pdata in set(plotData).intersection([RunTypes.SDC_Coarse, RunTypes.SDC_Fine]):
        numPlots += plotFile(plotType, root + '/' + pdata + '.dump')
    return numPlots


# dictionary containing the function pointer for each plot type
plotDict = {PlotTypes.LastIterationError: plotLastIterationError,
            PlotTypes.Iter: plotIterations,
            PlotTypes.LastStepError: plotLastStepError,
            PlotTypes.LastStepResidual: plotLastStepResidual,
            PlotTypes.AllIterationErrors: plotAllIterations,
            PlotTypes.AllIterationResiduals: plotAllIterations,
            PlotTypes.Timings: plotTimings}


def plotFile(plotType, fname):
    print(fname)
    if plotType in [PlotTypes.LastStepError,
                    PlotTypes.LastStepResidual,
                    PlotTypes.LastIterationError,
                    PlotTypes.AllIterationErrors,
                    PlotTypes.AllIterationResiduals,
                    PlotTypes.Iter,
                    PlotTypes.Timings]:
        if plotType not in [PlotTypes.Iter, PlotTypes.Timings]:
            pyplot.yscale('log')
        if plotType in [PlotTypes.AllIterationResiduals,
                        PlotTypes.LastStepResidual]:
            pyplot.ylabel('Residuum')
        elif plotType == PlotTypes.Iter:
            pyplot.ylabel('Iterationen')
        elif plotType == PlotTypes.Timings:
            pyplot.ylabel('Zeit[s]')
        else:
            pyplot.ylabel('Fehler')
    if plotType in [PlotTypes.LastStepError,
                    PlotTypes.LastStepResidual]:
        pyplot.xlabel('Iterationen')
    elif plotType == PlotTypes.Timings:
        pyplot.xlabel('Prozessor')
    else:
        pyplot.xlabel('Zeitschritt')

    with open(fname, 'rb') as input:
        dumpObj = pickle.load(input)

    data = dumpObj
    runtype = dumpObj.runtype
    if plotType == PlotTypes.AllIterationErrors:
        data = data.result.errMap
    elif plotType == PlotTypes.AllIterationResiduals:
        data = data.result.resMap
    return plotDict[plotType](data, getStyle(runtype, False), labelDict[runtype] + r" $N_p=$%d" % dumpObj.nproc)


def getSpeedUps(fileNames, sdcTime):
    data = {}
    files = glob.glob(fileNames)
    if len(files) == 0:
        return
    for file in files:
        print(file)
        with open(file, 'rb') as input:
            dumpObj = pickle.load(input)
            timeMeasure = numpy.max(dumpObj.result.timing.totalTimes)
        speedUp = sdcTime / timeMeasure
        efficiency = speedUp / dumpObj.nproc
        data[dumpObj.nproc] = [speedUp, efficiency]
        sortedData = sorted(data.items())
        nprocs = []
        speedUps = []
        efficiencies = []
        for t in sortedData:
            nprocs += [t[0]]
            speedUps += [t[1][0]]
            efficiencies += [t[1][1]]
        print(dumpObj.runtype, 'nproc:', dumpObj.nproc,
              'maxIter:', dumpObj.result.maxIter,
              'speedUp:', speedUp, 'efficiency:', efficiency)
    return nprocs, speedUps, efficiencies


def plotSpeedUp(plotEfficiency, plotData, root='.'):
    numPlots = 0
    with open('%s/%s.dump' % (root, RunTypes.SDC_Fine), 'rb') as input:
        dumpObj = pickle.load(input)
        sdcTime = dumpObj.result.timing.totalTimes[0]

    for type in list(speedUpSet & set(plotData)):
        ret = getSpeedUps('%s/%s_nproc*.dump' % (root, type), sdcTime)
        if ret is None:
            continue
        nprocs, speedUps, efficiency = ret
        if len(nprocs) > 0:
            if plotEfficiency:
                plotData = efficiency
            else:
                plotData = speedUps
            numPlots += 1
            pyplot.plot(nprocs, plotData, getStyle(type, False), label=labelDict[type])

    pyplot.xlabel(r'$N_p$')
    if plotEfficiency:
        pyplot.ylabel('Effizienz')
    else:
        pyplot.ylabel('Speedup')
    return numPlots


# this function plots the calculated theoretical Speedups for the juqueen runs, provided the file 'Timings' is located
# in the same folder. This file contains the timings of one step of the converged fine and coarse SDC propagator
def plotSpeedUpAlpha(plotData, root):
    numPlots = 0
    with open('%s/%s.dump' % (root, RunTypes.SDC_Fine), 'rb') as input:
        m = pickle.load(input).result.maxIter
        print(root, 'M=', m)
    with open('%s/%s.dump' % (root, RunTypes.SDC_Coarse), 'rb') as input:
        mc = pickle.load(input).result.maxIter
        print(root, 'Mc=', mc)
    with open('Timings', 'rb') as f:
        alphas = pickle.load(f)
    i = int(re.findall("runs_00(\d*)", root)[0])
    alphaClassic = alphas[i][0]/alphas[i][1]
    alphaHybrid = (alphas[i][0]/(mc+1))/(alphas[i][1]/(m+1))

    print('alphaClassic:', alphaClassic, ' alphaHybrid:', alphaHybrid)

    speedUpDict = {RunTypes.PARA_HYBRID_FULL: {},
                   RunTypes.PARA_HYBRID_PARTIAL: {},
                   RunTypes.PARA_CLASSIC: {}}
    fileNames = processPlotDataName(list(speedUpAlphaSet & set(plotData)), 0, root)
    for fname in fileNames:
        for file in glob.glob(fname):
            with open(file, 'rb') as input:
                dumpObj = pickle.load(input)
                nproc = dumpObj.nproc
                runType = dumpObj.runtype
                maxIter = dumpObj.result.maxIter
            if runType == RunTypes.PARA_CLASSIC:
                alpha = alphaClassic
                d = 1
            else:
                alpha = alphaHybrid
                d = m
            speedUp = d / (alpha + maxIter / nproc * (alpha + 1))
            speedUpDict[runType][nproc] = speedUp
    for k in speedUpDict.keys():
        numPlots += 1
        plotDictData(speedUpDict[k], labelDict[k] + ' Theorie', getStyle(k, True))
    return numPlots


def plotDictData(dict, label, style):
    sortedData = sorted(dict.items())
    np = []
    su = []
    for i in range(len(sortedData)):
        np += [sortedData[i][0]]
        su += [sortedData[i][1]]
    pyplot.plot(np, su, style, label=label)


def plotConvergenceData(plotData, root='.'):
    numPlots = 0
    fineOrders = []
    colorDict = {2: 'b', 3: 'g', 5: 'r'}
    for file in glob.glob('%s/orderPlot_numNodes[53].pkl' % root):
        with open(file, 'rb') as f:
            [errDict, dtArr, input] = pickle.load(f)

        for pdata in plotData:
            if pdata not in errDict.keys():
                continue
            errList = errDict[pdata]
            num_nodes = input.num_nodes
            if pdata == RunTypes.SDC_Coarse: num_nodes = input.num_nodes_coarse
            numPlots += 1

            if pdata not in [RunTypes.SDC_Coarse, RunTypes.SDC_Fine]:
                label = '%s #Knoten: %d, %d' % (labelDict[pdata], num_nodes, input.num_nodes_coarse)
            else:
                label = '%s #Knoten: %d' % (labelDict[pdata], num_nodes)

            pyplot.plot(dtArr, errList, getStyle(pdata, False), label=label,
                        color=colorDict[num_nodes])

        errListSDC = errDict[RunTypes.SDC_Fine]
        orderFine = round(abs((log10(errListSDC[1]) - log10(errListSDC[0])) / (log10(dtArr[1]) - log10(dtArr[0]))))
        fineOrders += [orderFine]
        p = log10(errListSDC[0]) - orderFine * abs(log10(dtArr[-1]) - log10(dtArr[0]))
        numPlots += 1
        pyplot.plot([dtArr[0], dtArr[-1]],
                    [errListSDC[0], pow(10, p)], '--',
                    label='%d-te Ordnung' % orderFine, color = colorDict[int((orderFine+2)/2)])
        errListSDC = errDict[RunTypes.SDC_Coarse]
        orderCoarse = round(abs((log10(errListSDC[1]) - log10(errListSDC[0])) / (log10(dtArr[1]) - log10(dtArr[0]))))
        if orderCoarse not in fineOrders:
            numPlots += 1
            p = log10(errListSDC[0]) - orderCoarse * abs(log10(dtArr[-1]) - log10(dtArr[0]))
            pyplot.plot([dtArr[0], dtArr[-1]],
                        [errListSDC[0], pow(10, p)], '--',
                        label='%d-te Ordnung' % orderCoarse, color = colorDict[int((orderCoarse+2)/2)])
        pyplot.xlim([dtArr[-1], dtArr[0]])
    pyplot.xscale('log')
    pyplot.yscale('log')
    ax = pyplot.gca()
    ax.invert_xaxis()

    pyplot.xlabel(r'$\Delta t$')
    pyplot.ylabel('Fehler')
    return numPlots


# main function of this script. It calls the corresponding plot function for each plot type in the list plotTypes.
# - The data of the runtypes in the list plotData will be plottet if available via the dumpfiles.
# - nprocs specifies for which number of processors the data will be plottet. If no value is provided the data will be
#   plottet for all available procs.
# - root defines the folder in which the dump files containing the data to be plottet are located.
def plot(plotTypes, plotData, nprocs=[0], root='.'):
    for pType in plotTypes:
        fig = pyplot.figure()
        ax = pyplot.subplot(111)

        numPlots = 0
        if pType in [PlotTypes.SpeedUp, PlotTypes.Efficiency, PlotTypes.SpeedUpAlpha]:
            if pType in [PlotTypes.SpeedUpAlpha]:
                numPlots += plotSpeedUpAlpha(plotData, root)
            numPlots += plotSpeedUp(pType == PlotTypes.Efficiency, plotData, root)
        elif pType == PlotTypes.ConvergenceData:
            numPlots = plotConvergenceData(plotData, root)
        else:
            numPlots = plotComparison(pType, plotData, nprocs, root)

        size = fig.get_size_inches() * fig.dpi
        factor = (math.ceil(numPlots/2) * 27)/size[1]
        box = ax.get_position()
        x0 = box.x0-box.width*0.02
        y0 = box.y0 + factor
        ax.set_position([x0, y0,
                         0.97-x0, 0.97-y0])

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.47, -box.y0 - box.height*0.02), ncol=2, prop={'size': 12})
        # pyplot.legend(loc=2)
    pyplot.show()
