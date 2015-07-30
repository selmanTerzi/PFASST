__author__ = 's.terzi'

import ProcessStarter
import pickle


class Input:
    dt = 0.01
    num_steps = 80
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
    PFASST = 'pfasst'

class DumpObj:
    def __init__(self, input, output, nproc):
        self.input = input
        self.output = output
        self.nproc = nproc

#----------------------------------------------------------------------------------------------------------------------#
#                                                   SDC                                                                #
#----------------------------------------------------------------------------------------------------------------------#
def run_vanilla_sdc(input):
    return ProcessStarter.run_vanilla_sdc(input.num_steps, input.dt, input.spatial_dofs, input.abs_res_tol,
                                          input.num_nodes, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#                                           PARAREAL_SERIAL                                                            #
#----------------------------------------------------------------------------------------------------------------------#
def run_parareal_serial(input):
    return ProcessStarter.run_parareal_classic(1, input.num_steps, input.dt, input.spatial_dofs,
                                               input.spatial_dofs_coarse, input.abs_res_tol, input.num_nodes,
                                               input.num_nodes_coarse, input.num_crse_iter, input.num_fine_iter)
#----------------------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------------------#
#                                            PARAREAL_MPI                                                              #
#----------------------------------------------------------------------------------------------------------------------#
def run_parareal_mpi(input, nproc):
    return ProcessStarter.run_parareal_classic(nproc, input.num_steps, input.dt, input.spatial_dofs,
                                               input.spatial_dofs_coarse, input.abs_res_tol,
                                               input.num_nodes, input.num_nodes_coarse, input.num_crse_iter,
                                               input.num_fine_iter, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------------------#
#                                           PARAREAL_HYBRID                                                            #
#----------------------------------------------------------------------------------------------------------------------#
def run_parareal_hybrid_full(input, nproc):
    return ProcessStarter.run_parareal_hybrid_full(nproc, input.num_steps, input.dt, input.spatial_dofs,
                                                   input.spatial_dofs_coarse, input.abs_res_tol, input.num_nodes,
                                                   input.num_nodes_coarse, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#                                                PFASST                                                                #
#----------------------------------------------------------------------------------------------------------------------#
def run_pfasst(input, nproc):
    return ProcessStarter.run_pfasst(nproc, input.num_steps, input.dt, input.spatial_dofs, input.spatial_dofs_coarse,
                                     input.abs_res_tol, input.num_nodes, input.num_nodes_coarse, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#

def dumpOutput(input, output, fileName, nproc = 1):
    with open(fileName, 'wb') as outfile:
        pickle.dump(DumpObj(input, output, nproc), outfile)