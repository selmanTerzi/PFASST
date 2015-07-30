__author__ = 's.terzi'

import ConvenientProcessStarter as cps
import pickle

input = cps.Input()

tend = 0.1
input.spatial_dofs = 128
input.spatial_dofs_coarse = 128
input.abs_res_tol = 1e-14
input.num_nodes = 5
input.num_nodes_coarse = 5
input.num_crse_iter = 1
input.num_fine_iter = 15

dtArr = []
errPara = []
resPara = []

input.num_nodes_coarse = input.num_nodes
i = 0

for input.dt in [1e-2, 1e-3]:
    input.num_steps = round(tend/input.dt)
    print('dt: %f num_steps %d' % (input.dt, input.num_steps))

    output = cps.run_parareal_serial(input)
    maxIter = output.maxIter
    maxStep = output.maxStep
    errIterList = []
    resIterList = []
    for i in range(1, maxIter+1):
        errIterList += [output.errMap[i][maxStep]]
        resIterList += [output.resMap[i][maxStep]]
    errPara += [errIterList]
    resPara += [resIterList]
    dtArr += [input.dt]

with open('convergencePlotData.pkl', 'wb') as output:
        pickle.dump([errPara, resPara, dtArr, input], output)