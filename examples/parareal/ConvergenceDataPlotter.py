__author__ = 's.terzi'

from matplotlib import pyplot
import glob
import pickle
from numpy import *

pyplot.figure()
legend = []
char = ['*', 'x', 's']
i = 0

for file in glob.glob('orderPlot_numNodes[1-9].pkl'):
    print(file)
    with open(file, 'rb') as input:
        [errPara, errSDC, dtArr, inputData] = pickle.load(input)

    order = abs((log10(errPara[1])-log10(errPara[0]))/(log10(dtArr[1])-log10(dtArr[0])))

    pyplot.plot(dtArr, errPara, '-s')
    pyplot.plot(dtArr, errSDC, '-o')

    realOrder = 2 * inputData.num_nodes - 2
    orderData = [10**(log10(dtArr[0])-1.2), 10**(log10(errPara[0])-realOrder*1.2)]
    pyplot.plot([dtArr[0], orderData[0]], [errPara[0], orderData[1]], '-.'+char[i])
    i = (i+1) % len(char)

    legend += ['parr, %d Knoten' % inputData.num_nodes,
               'sdc, %d Knoten' % inputData.num_nodes,
               '%dte Ordnung' % realOrder]

pyplot.xscale('log')
pyplot.yscale('log')
# pyplot.axis('equal')
ax = pyplot.gca()
ax.invert_xaxis()

pyplot.xlabel('dt')
pyplot.ylabel('Fehler')
pyplot.legend(legend)

pyplot.title('tend: %g spatial_dofs: %d abs_res_tol: %g'
             % (inputData.dt*inputData.num_steps, inputData.spatial_dofs, inputData.abs_res_tol))
pyplot.savefig('orderPlot.png')
pyplot.show()