__author__ = 's.terzi'

from matplotlib import pyplot
import pickle
import glob

def plotError(output, style, label):
    errMap = output['errMap']
    print(errMap)
    pyplot.plot(list(errMap[output['maxIter']].values()), style, label=label)

def plotIter(output, style, label):
    iterMap = output['iterMap']
    pyplot.plot(list(iterMap.values()), style, label=label)

# pyplot.figure()
# pyplot.yscale('log')


for file in glob.glob('*.dump'):
    print(file)
    with open(file, 'rb') as input:
        dumpObj = pickle.load(input)
    # plotError(dumpObj['output'], 's-', 'err_%s' % file)
    plotIter(dumpObj['output'], '^-', 'iter_%s' % file)

pyplot.legend()
pyplot.show()