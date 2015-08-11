__author__ = 's.terzi'
import re
import glob
from subprocess import *
from ProcessStarter import RunTypes

sciRegex = "\d*(?:\.\d+(?:[eE][-+]\d+)?)?"
timeString = "time Measurement"
errParaString = "ErrorParareal"

class Result:
    def __init__(self, errMap, resMap, iterMap, maxIter, maxStep, timeMeasure):
        self.errMap = errMap
        self.resMap = resMap
        self.iterMap = iterMap
        self.maxIter = maxIter
        self.maxStep = maxStep
        self.timeMeasure = timeMeasure

def getResultObj(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure):
    errMap, resMap = fillErrMapAndResMap(errMap, resMap, iterMap, maxIter, maxStep)
    return Result(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure)


def fillErrMapAndResMap(errMap, resMap, iterMap, maxIter, maxStep):
    for iter in range(1, maxIter+1):
        for step in range(1, maxStep+1):
            if step not in errMap[iter]:
                errMap[iter][step] = errMap[iter-1][step]
            if iter and iter - 1 in resMap and step not in resMap[iter]:
                resMap[iter][step] = resMap[iter-1][step]

    return errMap, resMap


def getLinesPFASST(params):
    file = params[0]
    spatial_dofs = params[1]
    grep = Popen(['grep', '\(n2:\s*%d\|%s\)' % (spatial_dofs, timeString)], stdin=file, stdout=PIPE)
    return grep.communicate()[0].splitlines()


def getLinesParaClassic(file):
    grep = Popen(['grep', '\(%s\|%s\)' % (errParaString, timeString)], stdin=file, stdout=PIPE)
    return grep.communicate()[0].splitlines()


def getLinesParaHybrid(file):
    grep2 = Popen(['grep', '\(step:\|%s\)' % timeString], stdin=PIPE, stdout=PIPE)
    grep1 = Popen(['grep', '-A1', '\(Fine Sweep\|%s\)' % timeString], stdin = file, stdout=grep2.stdin)
    output = grep2.communicate()[0]
    grep1.wait()
    return output.splitlines()


def getLinesSDC(file):
    grep = Popen(['grep', '\(step:\|%s\)' % timeString], stdin=file, stdout=PIPE)
    return grep.communicate()[0].splitlines()


def parseLines(lines):
    errMap = {}
    resMap = {}
    iterMap = {}

    maxStep = 0
    maxIter = 0
    timeMeasure = 0

    for i in range(len(lines)):
        line = lines[i].decode("utf-8")
        if len(re.findall(timeString, line)) > 0:
            print(line)
            t = float(re.findall("%s:\s*(%s)" % (timeString, sciRegex), line)[0])
            if t > timeMeasure: timeMeasure = t
        else:
            iter = int(re.findall("iter:\s*(\d*)", line)[0])
            step = int(re.findall("step:\s*(\d*)", line)[0])
            iterMap[step] = iter

            if step > maxStep: maxStep = step
            if iter > maxIter: maxIter = iter

            error = float(re.findall("err:\s*(%s)" % sciRegex, line)[0])
            errMap[iter] = errMap.get(iter, {})
            errMap[iter][step] = error

            residual = float(re.findall("residual:\s*(%s)" % sciRegex, line)[0])
            resMap[iter] = resMap.get(iter, {})
            resMap[iter][step] = residual

            # print("iter: %d step: %d error: %e" % (iter, step, error))

    print('timeMeasure: %g' % timeMeasure)
    return errMap, resMap, iterMap, maxIter, maxStep, timeMeasure

grepDict = {RunTypes.PARA_HYBRID_FULL: getLinesParaHybrid,
            RunTypes.PARA_HYBRID_PARTIAL: getLinesParaHybrid,
            RunTypes.PARA_CLASSIC: getLinesParaClassic,
            RunTypes.PARA_CLASSIC_SERIAL: getLinesParaClassic,
            RunTypes.PFASST: getLinesPFASST,
            RunTypes.SDC_Coarse: getLinesSDC,
            RunTypes.SDC_Fine: getLinesSDC}

def printMap(map):
    for k in sorted(map):
        for n in sorted(map[k]):
            print("iter: %d step: %d error: %e" % (k, n, map[k][n]))

def getResult(runtype, input, dir='.'):

    errMap = {}
    resMap = {}
    iterMap = {}
    maxIter = 0
    maxStep = 0
    timeMeasure = 0

    pattern = dir + '/*.log'
    for fileName in sorted(glob.glob(pattern)):
        print(fileName)
        with open(fileName, 'r') as file:
            param = file
            if runtype == RunTypes.PFASST:
                param = param, input.spatial_dofs
            errMapLocal, resMapLocal, iterMapLocal, maxIterLocal, maxStepLocal, t = parseLines(grepDict[runtype](param))

        mergeDictionaries(errMap, errMapLocal)
        mergeDictionaries(resMap, resMapLocal)
        iterMap.update(iterMapLocal)

        if maxIterLocal > maxIter: maxIter = maxIterLocal
        if maxStepLocal > maxStep: maxStep = maxStepLocal
        if t > timeMeasure: timeMeasure = t

    return getResultObj(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure)


def mergeDictionaries(dicA, dicB):
    for k in sorted(dicB):
        if k not in dicA.keys():
            dicA[k] = dicB[k]
        else:
            dicA[k].update(dicB[k])