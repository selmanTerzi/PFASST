__author__ = 's.terzi'

import os
import subprocess
from ProcessStarter import Input
from OutputDumper import *
import pickle

input = Input()
input.spatial_dofs = 4096
input.num_nodes = 5
input.abs_res_tol = 1e-10
input.diffRes = True

subprocess.call("rm -rf *.dump", shell=True)
runTypes = [
            RunTypes.SDC_Fine,
            RunTypes.SDC_Coarse,
            # RunTypes.PFASST,
            # RunTypes.PARA_HYBRID_FULL,
            # RunTypes.PARA_HYBRID_PARTIAL,
            # RunTypes.PARA_CLASSIC,
            ]


def getRuntimes():
    with open('sdc_coarse.dump', 'rb') as f:
        tc = pickle.load(f).result.timeMeasure
    with open('sdc_fine.dump', 'rb') as f:
        tf = pickle.load(f).result.timeMeasure
    return tc, tf

def getAlpha(input):
    subprocess.call("rm -rf *.dump", shell=True)
    input.num_iter = 1
    dumpRuns(runTypes, 0, input)
    tmCoarse, tmFine = getRuntimes()

    subprocess.call("rm -rf *.dump", shell=True)
    input.num_iter = 2
    dumpRuns(runTypes, 0, input)
    tc, tf = getRuntimes()
    alphaHybrid = (tc - tmCoarse) / (tf - tmFine)

    subprocess.call("rm -rf *.dump", shell=True)
    input.num_iter = 20
    dumpRuns(runTypes, 0, input)
    tc, tf = getRuntimes()
    alphaClassic = tc / tf
    return alphaHybrid, alphaClassic

def getAlphas():
    input.tend = 0.005
    input.dt = 0.005
    input.num_nodes_coarse = (input.num_nodes + 1)/2
    input.spatial_dofs_coarse = input.spatial_dofs/2

    alpha1Hybrid, alpha1Classic = getAlpha(input)

    input.tend = 0.01
    input.dt = 0.01

    alpha2Hybrid, alpha2Classic = getAlpha(input)

    input.num_nodes_coarse = input.num_nodes

    alpha3Hybrid, alpha3Classic = getAlpha(input)

    input.num_nodes_coarse = (input.num_nodes + 1)/2
    input.spatial_dofs_coarse = input.spatial_dofs

    alpha4Hybrid, alpha4Classic = getAlpha(input)

    outObj = [alpha1Classic, alpha1Hybrid,
              alpha2Classic, alpha2Hybrid,
              alpha3Classic, alpha3Hybrid,
              alpha4Classic, alpha4Hybrid]

    return outObj

rep = 200
outObj = [0, 0, 0, 0, 0, 0, 0, 0]
for i in range(rep):
    alphas = getAlphas()
    for j in range(len(alphas)):
        outObj[j] += alphas[j]


for j in range(len(alphas)):
    outObj[j] /= rep
print(outObj)

with open('Alphas.al', 'wb') as f:
    pickle.dump(outObj, f)
