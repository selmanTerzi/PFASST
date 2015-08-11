__author__ = 's.terzi'

import os
import subprocess
from ProcessStarter import Input
from OutputDumper import *

input = Input()
input.spatial_dofs = 4096
input.spatial_dofs_coarse = input.spatial_dofs/2
input.num_nodes = 5
input.num_nodes_coarse = (input.num_nodes + 1)/2
input.abs_res_tol = 1e-10
input.tend = 0.64
input.dt = 0.01

# subprocess.call("rm -rf *.dump", shell=True)
# runTypes = [
#             RunTypes.SDC_Fine,
#             # RunTypes.SDC_Coarse,
#             # RunTypes.PFASST,
#             # RunTypes.PARA_HYBRID_FULL,
#             # RunTypes.PARA_HYBRID_PARTIAL,
#             # RunTypes.PARA_CLASSIC,
#             ]
# nprocs = [
#           2,
#           4,
#           8
#           ]
#
# dumpRuns(runTypes, nprocs, input)

dirs = ['paraClassic*', 'paraPartial*', 'pfasst_*']
for dir in dirs:
    dumpDir(dir)