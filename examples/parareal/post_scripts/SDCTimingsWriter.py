import subprocess
import pickle
import re
from OutputDumper import *

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


def getTimings(roots):
    timings = {}
    for root in roots:
        dir = runDirectorysRoot + '/' + root + '/alphas'
        i = int(re.findall("runs_00(\d*)", root)[0])
        tCoarse = getRuntime(dir + '/sdc_coarse/*.log')
        tFine = getRuntime(dir + '/sdc_fine/*.log')
        timings[i] = [tCoarse, tFine]
    return timings


roots = [
        'runs_002',
        'runs_003',
        'runs_004'
        ]
with open('Timings', 'wb') as f:
    pickle.dump(getTimings(roots), f)
