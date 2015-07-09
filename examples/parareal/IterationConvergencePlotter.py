__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
from numpy import *

pyplot.figure()
legend = []
char = ['*', 'x', 's']
i = 0

with open('convergencePlotData.pkl', 'rb') as input:
    [errPara, resPara, dtArr, inputData] = pickle.load(input)

for i in range(len(errPara)):
    iterations = arange(1, len(errPara[i])+1, 1)
    pyplot.plot(iterations, errPara[i], '-s')
    pyplot.plot(iterations, resPara[i], '-.o')

    legend += ['err dt: %g' % dtArr[i],
               'res dt: %g' % dtArr[i]]

pyplot.yscale('log')

pyplot.xlabel('Parareal-Iteration')
pyplot.legend(legend)

pyplot.title('tend: %g spatial_dofs: %d num_nodes: %d abs_res_tol: %g'
             % (inputData.dt * inputData.num_steps, inputData.spatial_dofs, inputData.num_nodes, inputData.abs_res_tol))
pyplot.savefig('convergencePlot.png')
pyplot.show()