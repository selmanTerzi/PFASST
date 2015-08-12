__author__ = 's.terzi'

import ProcessStarter as ps
from ProcessStarter import *
from OutputParser import getResult
import os
import glob
import pickle
from subprocess import *
import re
import copy

runDirectorysRoot = os.path.expanduser("~")+'/project/PFASST/run/'

class DumpObj:
    def __init__(self, result, nproc):
        self.result = result
        self.nproc = nproc


def dumpOutput(output, fileName, nproc = 1):
    with open(fileName, 'wb') as outfile:
        pickle.dump(DumpObj(input, output, nproc), outfile)


def getCoarseInput(input):
    inputCoarse = copy.deepcopy(input)
    inputCoarse.spatial_dofs = input.spatial_dofs_coarse
    inputCoarse.num_nodes = input.num_nodes_coarse
    return inputCoarse


def dumpRuns(runTypes, nprocs, inputFine):
    for runType in runTypes:
        if runType in [RunTypes.SDC_Fine, RunTypes.SDC_Coarse]:
            if runType == RunTypes.SDC_Coarse:
                inputCurrent = getCoarseInput(inputFine)
            else:
                inputCurrent = inputFine
            ps.run_prog(runType, inputCurrent)
            dumpOutput(getResult(runType, inputCurrent), '%s.dump' % runType)
        else:
            for nproc in nprocs:
                ps.run_prog(runType, inputFine, nproc)
                dumpOutput(getResult(runType, inputFine), '%s_nproc%d.dump' % (runType, nproc), nproc)


def getRunTypeAndInput(dir):
    with open(glob.glob(dir + '/*.job')[0]) as f:
        grep = Popen(['grep', '\(--np\|/PFASST/\|--tend\)'], stdin=f, stdout=PIPE)
        output = grep.communicate()[0].splitlines()

    procName = []
    np = []
    tend = []
    dt = []
    sdofs = []
    sdofsC = []
    nodes = []
    nodesC = []

    for i in range(len(output)):
        line = output[i].decode("utf-8")
        procName += re.findall('/([^/]*)\s+\\\\', line)
        np += re.findall('--np\s*(\d*)', line)
        tend += re.findall('--tend\s*(\d*)', line)
        dt += re.findall('--dt\s*(\d*)', line)
        sdofs += re.findall('--spatial_dofs\s*(\d*)', line)
        sdofsC += re.findall('--spatial_dofs_coarse\s*(\d*)', line)
        nodes += re.findall('--num_nodes\s*(\d*)', line)
        nodesC += re.findall('--num_nodes_coarse\s*(\d*)', line)

    runType = progName2RunType[procName[0]]
    input = ps.Input()
    nproc = int(np[0])
    input.tend = float(tend[0])
    input.dt = float(dt[0])
    input.spatial_dofs = int(sdofs[0])
    input.num_nodes = int(nodes[0])
    if len(sdofsC) > 0: input.spatial_dofs_coarse = int(sdofsC[0])
    if len(nodesC) > 0: input.num_nodes_coarse = int(nodesC[0])

    return runType, input, nproc


def dumpDir(dirNames):
    for dir in glob.glob(runDirectorysRoot + dirNames):
        print(dir)
        runType, input, nproc = getRunTypeAndInput(dir)
        if runType == RunTypes.SDC_Fine:
            fname = '%s.dump' % RunTypes.SDC_Fine
        else:
            fname = '%s_nproc%d.dump' % (runType, nproc)
        print(fname)
        dumpOutput(getResult(runType, input, dir), fname, nproc)

