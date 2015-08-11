__author__ = 's.terzi'

import os
from subprocess import *

buildFolder = os.path.expanduser("~")+'/project/PFASST/build_mpi/'
examplesParareal = buildFolder + 'examples/parareal/'
examplesAdvecDiffusion = buildFolder + 'examples/advection_diffusion/'


class Input:
    dt = 0.01
    tend = 0.1
    spatial_dofs = 128
    spatial_dofs_coarse = 64
    num_nodes = 5
    num_nodes_coarse = 3
    abs_res_tol = 1e-10
    num_crse_iter = 20
    num_fine_iter = 20
    num_iter = 20


class RunTypes:
    SDC_Fine = 'sdc_fine'
    SDC_Coarse = 'sdc_coarse'
    PARA_CLASSIC_SERIAL = 'parareal_classic_serial'
    PARA_CLASSIC = 'parareal_classic'
    PARA_HYBRID_FULL = 'parareal_hybrid_full'
    PARA_HYBRID_PARTIAL = 'parareal_hybrid_partial'
    PFASST = 'pfasst'

progName2RunType = {'vanilla_sdc': RunTypes.SDC_Fine,
                    'parareal_classic': RunTypes.PARA_CLASSIC,
                    'parareal_hybrid_full': RunTypes.PARA_HYBRID_FULL,
                    'parareal_hybrid_partial': RunTypes.PARA_HYBRID_PARTIAL,
                    'mpi_pfasst': RunTypes.PFASST}

progDict = {RunTypes.SDC_Coarse: examplesAdvecDiffusion+'vanilla_sdc',
            RunTypes.SDC_Fine: examplesAdvecDiffusion+'vanilla_sdc',
            RunTypes.PARA_CLASSIC_SERIAL: examplesParareal+RunTypes.PARA_CLASSIC_SERIAL,
            RunTypes.PARA_CLASSIC: examplesParareal+RunTypes.PARA_CLASSIC,
            RunTypes.PARA_HYBRID_FULL: examplesParareal+RunTypes.PARA_HYBRID_FULL,
            RunTypes.PARA_HYBRID_PARTIAL: examplesParareal+RunTypes.PARA_HYBRID_PARTIAL,
            RunTypes.PFASST: examplesAdvecDiffusion+'mpi_pfasst'}


def run_prog(runType, input, nproc=1):
    prog = progDict[runType]
    pargs = []
    if nproc > 1:
        pargs = ['mpirun', '-np', '%d' % nproc]

    pargs += [prog, '-q', '-c',
              '--dt', '%f' % input.dt,
              '--tend', '%f' % input.tend,
              '--spatial_dofs', '%d' % input.spatial_dofs,
              '--num_nodes', '%d' % input.num_nodes,
              '--num_iter', '%d' % input.num_iter,
              '--abs_res_tol', '%g' % input.abs_res_tol]
    if runType not in [RunTypes.SDC_Coarse, RunTypes.SDC_Fine]:
        pargs += ['--spatial_dofs_coarse', '%d' % input.spatial_dofs_coarse,
                 '--num_nodes_coarse', '%d' % input.num_nodes_coarse]
    if runType in [RunTypes.PARA_CLASSIC, RunTypes.PARA_CLASSIC_SERIAL]:
        pargs += ['--num_crse_iter', '%d' % input.num_crse_iter,
                 '--num_fine_iter', '%d' % input.num_fine_iter]

    call("rm -rf *.log", shell=True)
    with open(os.devnull, "w") as f:
        call(pargs, stdout=f)
    print('%s_nproc%s' % (runType, nproc))