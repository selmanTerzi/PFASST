__author__ = 's.terzi'

import ConvenientProcessStarter as cps
from ConvenientProcessStarter import RunTypes
import subprocess

input = cps.Input()
input.spatial_dofs = 64
input.spatial_dofs_coarse = input.spatial_dofs/2
input.num_nodes = 5
input.num_nodes_coarse = (input.num_nodes + 1)/2
input.num_iter = 10
input.num_steps = 40
input.abs_res_tol = 1e-10

def dumpRuns(runTypes, nprocs):
    if RunTypes.SDC_Fine in runTypes:
        cps.dumpOutput(input, cps.run_vanilla_sdc(input), '%s.dump' % RunTypes.SDC_Fine)
    for nproc in nprocs:
        if RunTypes.PARA_CLASSIC in runTypes:
            cps.dumpOutput(input, cps.run_parareal_mpi(input, nproc), '%s_nproc%d.dump' % (RunTypes.PARA_CLASSIC,
                                                                                           nproc), nproc)
        if RunTypes.PARA_HYBRID in runTypes:
            cps.dumpOutput(input, cps.run_parareal_hybrid(input, nproc), '%s_nproc%d.dump' % (RunTypes.PARA_HYBRID,
                                                                                              nproc), nproc)
        if RunTypes.PFASST in runTypes:
            cps.dumpOutput(input, cps.run_pfasst(input, nproc), '%s_nproc%d.dump' % (RunTypes.PFASST,
                                                                                     nproc), nproc)
    if RunTypes.SDC_Coarse in runTypes:
        input.num_nodes = input.num_nodes_coarse
        input.spatial_dofs = input.spatial_dofs_coarse
        cps.dumpOutput(input, cps.run_vanilla_sdc(input), '%s.dump' % RunTypes.SDC_Coarse)

subprocess.call("rm -rf *.dump", shell=True)

runTypes = [
            RunTypes.SDC_Fine,
            RunTypes.SDC_Coarse,
            RunTypes.PFASST,
            RunTypes.PARA_HYBRID,
            RunTypes.PARA_CLASSIC,
            ]
nprocs = [2, 4, 8]

dumpRuns(runTypes, nprocs)