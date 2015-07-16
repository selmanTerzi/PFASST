__author__ = 's.terzi'

import ConvenientProcessStarter as cps
import subprocess

input = cps.Input()
input.spatial_dofs = 128
input.spatial_dofs_coarse = 32
input.num_nodes_coarse = 3
input.num_iter = 50
input.num_steps = 40

def dumpRuns(RunTypes, nprocs):
    if 'SDC' in RunTypes:
        cps.dumpOutput(input, cps.run_vanilla_sdc(input), 'sdc.dump')
    for nproc in nprocs:
        if 'PARA' in RunTypes:
            cps.dumpOutput(input, cps.run_parareal_mpi(input, nproc), 'parareal_mpi_nproc%d.dump' % nproc, nproc)
        if 'PARA_HYBRID' in RunTypes:
            cps.dumpOutput(input, cps.run_parareal_hybrid(input, nproc), 'parareal_hybrid_nproc%d.dump' % nproc, nproc)
        if 'PFASST' in RunTypes:
            cps.dumpOutput(input, cps.run_pfasst(input, nproc), 'pfasst_nproc%d.dump' % nproc, nproc)

subprocess.call("rm -rf *.dump", shell=True)

runTypes = [
            'SDC',
            # 'PARA',
            'PARA_HYBRID',
            'PFASST'
            ]
nprocs = [2, 4, 5, 8]

dumpRuns(runTypes, nprocs)


