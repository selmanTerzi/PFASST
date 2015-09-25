__author__ = 's.terzi'

import ProcessStarter as ps
from ProcessStarter import RunTypes
from OutputParser import getResult
import pickle

input = ps.Input()
tend = 0.2
input.spatial_dofs = 128
input.spatial_dofs_coarse = 64
input.abs_res_tol = 1e-14


def getLastError(runType, nproc=1):
    ps.run_prog(runType, input, nproc)
    result = getResult(runType, input)
    return result.errMap[result.maxIter][result.maxStep]

def runConvergenceTest(runTypes, num_nodes, dtSteps):
    for node in num_nodes:
        dtArr = []
        errDict = {}

        input.num_nodes = node
        input.num_nodes_coarse = node
        i = 0

        for i in range(dtSteps):
            input.dt = tend/2**(i+1)
            print('num_nodes: %d dt: %f' % (input.num_nodes, input.dt))

            for runType in runTypes:
                errList = errDict.get(runType, [])
                if runType not in [RunTypes.SDC_Fine, RunTypes.PARA_CLASSIC_SERIAL]:
                    nproc = 4
                else:
                    nproc = 1
                errList += [getLastError(runType, nproc)]
                errDict[runType] = errList

            dtArr += [input.dt]

        with open('orderPlot_numNodes%d.pkl' % node, 'wb') as output:
            pickle.dump([errDict, dtArr, input], output)

runTypes = [
            RunTypes.SDC_Fine,
            RunTypes.SDC_Coarse,
            RunTypes.PARA_CLASSIC,
            RunTypes.PARA_HYBRID_FULL,
            RunTypes.PARA_HYBRID_PARTIAL
            ]

runConvergenceTest(runTypes, [3, 5], 5)