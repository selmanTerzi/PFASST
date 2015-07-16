__author__ = 's.terzi'

import ConvenientProcessStarter as cps
import pickle

input = cps.Input()

tend = 0.2
input.spatial_dofs = 128
input.spatial_dofs_coarse = 64
input.abs_res_tol = 1e-14
input.num_nodes = 3
input.num_nodes_coarse = 3
input.num_crse_iter = 1
input.num_fine_iter = 15


for num_nodes in [3, 5]:

    dtArr = []
    errPara = []
    errSDC = []

    input.num_nodes = num_nodes
    input.num_nodes_coarse = num_nodes
    i = 0

    for i in range(6):
        input.dt = tend/2**(i+1)
        input.num_steps = round(tend/input.dt)
        print('num_nodes: %d dt: %f num_steps: %d' % (input.num_nodes, input.dt, input.num_steps))

        output = cps.run_parareal_serial(input)

        errPara += [output.errMap[output.maxIter][output.maxStep]]

        output = cps.run_vanilla_sdc(input)

        errSDC += [output.errMap[output.maxIter][output.maxStep]]

        dtArr += [input.dt]

    with open('orderPlot_numNodes%d.pkl' % num_nodes, 'wb') as output:
        pickle.dump([errPara, errSDC, dtArr, input], output)