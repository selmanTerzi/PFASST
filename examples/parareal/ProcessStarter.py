__author__ = 's.terzi'

import re
import glob
import os
from numpy import *
from subprocess import *

sciRegex = "\d*(?:\.\d+(?:[eE][-+]\d+)?)?"
buildFolder = "/home0/s.terzi/project/PFASST/build_mpi/"
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

def run_parareal_classic(nproc=2,
                         step=5,
                         dt=0.01,
                         spatial_dofs=4096,
                         spatial_dofs_coarse=4096,
                         abs_res_tol=1e-12,
                         num_nodes=3,
                         num_nodes_coarse=3,
                         num_crse_iter=20,
                         num_fine_iter=20,
                         num_iter=20):
    prog = 'parareal_classic'
    pargs = []
    if nproc == 1:
        prog += '_serial'
    else:
        pargs = ['mpirun', '-np', '%d' % nproc]


    grep = Popen(['grep', '\(%s\|%s\)' % (errParaString, timeString)], stdin=PIPE, stdout=PIPE)
    para = Popen(pargs + ['%s/examples/parareal/%s' % (buildFolder, prog),
                  '--dt', '%f' % dt,
                  '--tend', '%f' % (dt * step),
                  '--spatial_dofs', '%d' % spatial_dofs,
                  '--spatial_dofs_coarse', '%d' % spatial_dofs_coarse,
                  '--num_nodes', '%d' % num_nodes,
                  '--num_nodes_coarse', '%d' % num_nodes_coarse,
                  '--abs_res_tol', '%e' % abs_res_tol,
                  '--num_crse_iter', '%d' % num_crse_iter,
                  '--num_fine_iter', '%d' % num_fine_iter,
                  '--num_par_iter', '%d' % num_iter],
                 stdout=grep.stdin)
    output = grep.communicate()[0]
    para.wait()

    print('parareal_classic_nproc%s' % nproc)
    return getResultSDC(output)


def run_vanilla_sdc(step=5,
                    dt=0.01,
                    spatial_dofs=4096,
                    abs_res_tol=1e-12,
                    num_nodes=3,
                    num_iter=5):
    grep = Popen(['grep', '\(step:\|%s\)' % timeString], stdin=PIPE, stdout=PIPE)
    sdc = Popen(['%s/examples/advection_diffusion/vanilla_sdc' % buildFolder,
                 '--dt', '%f' % dt,
                 '--tend', '%f' % (dt * step),
                 '--spatial_dofs', '%d' % spatial_dofs,
                 '--abs_res_tol', '%e' % abs_res_tol,
                 '--num_iter', '%d' % num_iter,
                 '--num_nodes', '%d' % num_nodes],
                stdout=grep.stdin)
    output = grep.communicate()[0]
    sdc.wait()

    print('vanilla_sdc')
    return getResultSDC(output)


def run_pfasst(nproc=2,
               step=5,
               dt=0.01,
               spatial_dofs=4096,
               spatial_dofs_coarse=4096,
               abs_res_tol=1e-12,
               num_nodes=3,
               num_nodes_coarse=3,
               num_iter=20):
    call("rm -rf *.log", shell=True)
    grep = Popen(['grep', '\(n2:\s*%d\|%s\)' % (spatial_dofs, timeString)], stdin=PIPE, stdout=PIPE)
    pfasst = Popen(['mpirun', '-np', '%d' % nproc, '%s/examples/advection_diffusion/mpi_pfasst' % buildFolder,
                    '--dt', '%f' % dt,
                    '--tend', '%f' % (dt * step),
                    '--spatial_dofs', '%d' % spatial_dofs,
                    # '--spatial_dofs_coarse', '%d' % spatial_dofs_coarse,
                    '--num_nodes', '%d' % num_nodes,
                    # '--num_nodes_coarse', '%d' % num_nodes_coarse,
                    '--abs_res_tol', '%e' % abs_res_tol,
                    '--num_iter', '%d' % num_iter], stdout=grep.stdin)
    output = grep.communicate()[0]
    pfasst.wait()

    print('pfasst_nproc%s' % nproc)
    return getResultSDC(output)


def run_parareal_hybrid_full(nproc=2,
                            step=5,
                            dt=0.01,
                            spatial_dofs=4096,
                            spatial_dofs_coarse=4096,
                            abs_res_tol=1e-12,
                            num_nodes=3,
                            num_nodes_coarse=3,
                            num_iter=20):

    call("rm -rf *.log", shell=True)
    with open(os.devnull, "w") as f:
        call(['mpirun', '-np', '%d' % nproc, '%s/examples/parareal/parareal_hybrid' % buildFolder,
          '--dt', '%f' % dt,
          '--tend', '%f' % (dt * step),
          '--spatial_dofs', '%d' % spatial_dofs,
          '--spatial_dofs_coarse', '%d' % spatial_dofs_coarse,
          '--num_nodes', '%d' % num_nodes,
          '--num_nodes_coarse', '%d' % num_nodes_coarse,
          '--abs_res_tol', '%e' % abs_res_tol,
          '--num_iter', '%d' % num_iter], stdout=f)

    print('parareal_hybrid_full_nproc%s' % nproc)
    return getResultParaHybrid()


def getResultSDC(output):
    errMap, resMap, iterMap, maxIter, maxStep, timeMeasure = parseLines(output.splitlines())
    return getResult(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure)


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

    print('timeMeasure: %g' % timeMeasure)
    return errMap, resMap, iterMap, maxIter, maxStep, timeMeasure


def getResult(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure):
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


def getResultParaHybrid():

    errMap = {}
    resMap = {}
    iterMap = {}
    maxIter = 0
    maxStep = 0
    timeMeasure = 0

    for file in glob.glob('*.log'):
        grep2 = Popen(['grep', '\(step:\|%s\)'%timeString], stdin=PIPE, stdout=PIPE)
        grep1 = Popen(['grep', '-A1', '\(Fine Sweep\|%s\)' % timeString, '%s' % file], stdout=grep2.stdin)
        output = grep2.communicate()[0]
        grep1.wait()

        errMapLocal, resMapLocal, iterMapLocal, maxIterLocal, maxStepLocal, t = parseLines(output.splitlines())
        mergeDictionaries(errMap, errMapLocal)
        mergeDictionaries(resMap, resMapLocal)
        iterMap.update(iterMapLocal)

        if maxIterLocal > maxIter: maxIter = maxIterLocal
        if maxStepLocal > maxStep: maxStep = maxStepLocal
        if t > timeMeasure: timeMeasure = t

    return getResult(errMap, resMap, iterMap, maxIter, maxStep, timeMeasure)


def mergeDictionaries(dicA, dicB):
    for k in dicB.keys():
        if k not in dicA.keys():
            dicA[k] = dicB[k]
        else:
            dicA[k].update(dicB[k])