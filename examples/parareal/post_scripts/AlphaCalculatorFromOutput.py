__author__ = 's.terzi'

import subprocess
from OutputDumper import *
import pickle
import re

timeString = "time Measurement"
sciRegex = "\d*(?:\.\d+(?:[eE][-+]\d+)?)?"

def getRuntime(fileName):
    with open(glob.glob(fileName)[0], 'r') as file:
        grep = Popen(['grep', 'time Measurement'], stdin=file, stdout=PIPE)
    lines = grep.communicate()[0].splitlines()
    t = 0
    for line in lines:
        line = line.decode("utf-8")
        t += float(re.findall("%s:\s*(%s)" % (timeString, sciRegex), line)[0])
    print(fileName, len(lines))
    t /= len(lines)
    return t

def getAlphas(roots):
    alphas = {}
    for root in roots:
        dir = runDirectorysRoot + '/' + root + '/alphas'
        i = int(re.findall("runs_00(\d*)", root)[0])
        tCoarsePredict = getRuntime(dir + '/sdc_coarse_predict/*.log')
        tFinePredict = getRuntime(dir + '/sdc_fine_predict/*.log')
        tCoarseHybrid = getRuntime(dir + '/sdc_coarse_hybrid/*.log')
        tFineHybrid = getRuntime(dir + '/sdc_fine_hybrid/*.log')
        tCoarseClassic = getRuntime(dir + '/sdc_coarse_classic/*.log')
        tFineClassic = getRuntime(dir + '/sdc_fine_classic/*.log')
        alphas[i] = [tCoarseClassic/tFineClassic,
                     (tCoarseHybrid-tCoarsePredict)/(tFineHybrid-tFinePredict)]
    print(alphas)
    return alphas

roots = [
        'runs_002',
        'runs_003',
        'runs_004'
        ]
with open('Alphas.al', 'wb') as f:
    pickle.dump(getAlphas(roots), f)
