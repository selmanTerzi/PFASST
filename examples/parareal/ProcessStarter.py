__author__ = 's.terzi'

import re
import glob
import os
from numpy import *
from subprocess import *

sciRegex = "\d*(?:\.\d+(?:[eE][-+]\d+)?)?"
buildFolder = "/home0/s.terzi/project/PFASST/build_mpi/"


def run_parareal_serial(step=5,
                        dt=0.01,
                        spatial_dofs=4096,
                        spatial_dofs_coarse=4096,
                        abs_res_tol=1e-12,
                        num_nodes=3,
                        num_nodes_coarse=3,
                        num_crse_iter=20,
                        num_fine_iter=20,
                        num_iter=20):
    grep = Popen(['grep', 'ErrorParareal:'], stdin=PIPE, stdout=PIPE)
    para = Popen(['%s/examples/parareal/parareal_serial' % buildFolder,
                  '--dt', '%f' % dt,
                  '--num_steps', '%d' % step,
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

    return getErrorMapSDC(output)


def run_parareal_mpi(nproc=2,
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
    grep = Popen(['grep', 'ErrorParareal:'], stdin=PIPE, stdout=PIPE)
    para = Popen(['mpirun',
                  '-np', '%d' % nproc,
                  '%s/examples/parareal/parareal_mpi' % buildFolder,
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

    return getErrorMapSDC(output)

def run_vanilla_sdc(step=5,
                    dt=0.01,
                    spatial_dofs=4096,
                    abs_res_tol=1e-14,
                    num_nodes=3,
                    num_iter=5):
    grep = Popen(['grep', 'step:'], stdin=PIPE, stdout=PIPE)
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

    return getErrorMapSDC(output)

def getErrorMapSDC(output):
    lines = output.splitlines()
    errMap, resMap, iterMap, maxIter, maxStep = getErrAndResMap(lines)
    return getResult(errMap, resMap, iterMap, maxIter, maxStep)

def getErrAndResMap(lines):
    errMap = {}
    resMap = {}
    iterMap = {}

    maxStep = 0
    maxIter = 0

    for i in range(len(lines)):
        line = lines[i].decode("utf-8")
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

    return errMap, resMap, iterMap, maxIter, maxStep

def fillErrMapAndResMap(errMap, resMap, iterMap, maxIter, maxStep):
    for iter in range(1, maxIter+1):
        for step in range(1, maxStep+1):
            if step not in errMap[iter]:
                errMap[iter][step] = errMap[iter-1][step]
            if iter and iter - 1 in resMap and step not in resMap[iter]:
                resMap[iter][step] = resMap[iter-1][step]

    return errMap, resMap, iterMap

def getResult(errMap, resMap, iterMap, maxIter, maxStep):
    errMap, resMap, iterMap = fillErrMapAndResMap(errMap, resMap, iterMap, maxIter, maxStep)

    result = {}
    result['errMap'] = errMap
    result['resMap'] = resMap
    result['iterMap'] = iterMap
    result['maxIter'] = maxIter
    result['maxStep'] = maxStep

    return result

def run_parareal_hybrid(nproc=2,
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

    return getErrorMapParaHybrid()

def getErrorMapParaHybrid():

    errMap = {}
    resMap = {}
    iterMap = {}
    maxIter = 0
    maxStep = 0

    for file in glob.glob('*.log'):
        grep2 = Popen(['grep', 'step:'], stdin=PIPE, stdout=PIPE)
        grep1 = Popen(['grep', '-A1', 'Fine Sweep', '%s'%file], stdout=grep2.stdin)
        output = grep2.communicate()[0]
        grep1.wait()

        lines = output.splitlines()

        errMapLocal, resMapLocal, iterMapLocal, maxIterLocal, maxStepLocal = getErrAndResMap(lines)
        mergeDictionaries(errMap, errMapLocal)
        mergeDictionaries(resMap, resMapLocal)
        iterMap.update(iterMapLocal)

        if maxIterLocal > maxIter: maxIter = maxIterLocal
        if maxStepLocal > maxStep: maxStep = maxStepLocal


    return getResult(errMap, resMap, iterMap, maxIter, maxStep)

def mergeDictionaries(dicA, dicB):
    for k in dicB.keys():
        if k not in dicA.keys():
            dicA[k] = dicB[k]
        else:
            dicA[k].update(dicB[k])