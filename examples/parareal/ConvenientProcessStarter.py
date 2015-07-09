__author__ = 's.terzi'

import ProcessStarter
import pickle

class Input():
    dt = 0.01
    num_steps = 75
    spatial_dofs = 128
    spatial_dofs_coarse = 64
    num_nodes = 5
    num_nodes_coarse = 3
    abs_res_tol = 1e-14
    num_crse_iter = 20
    num_fine_iter = 20
    num_iter = 20

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
    return ProcessStarter.run_parareal_serial(input.num_steps, input.dt, input.spatial_dofs, input.spatial_dofs_coarse,
                                              input.abs_res_tol, input.num_nodes, input.num_nodes_coarse,
                                              input.num_crse_iter, input.num_fine_iter)
#----------------------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------------------#
#                                            PARAREAL_MPI                                                              #
#----------------------------------------------------------------------------------------------------------------------#
def run_parareal_mpi(input, nproc):
    return ProcessStarter.run_parareal_mpi(nproc, input.num_steps, input.dt, input.spatial_dofs,
                                           input.spatial_dofs_coarse, input.abs_res_tol,
                                           input.num_nodes, input.num_nodes_coarse, input.num_crse_iter,
                                           input.num_fine_iter, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------------------#
#                                           PARAREAL_HYBRID                                                            #
#----------------------------------------------------------------------------------------------------------------------#
def run_parareal_hybrid(input, nproc):
    return ProcessStarter.run_parareal_hybrid(nproc, input.num_steps, input.dt, input.spatial_dofs,
                                              input.spatial_dofs_coarse, input.abs_res_tol, input.num_nodes,
                                              input.num_nodes_coarse, input.num_iter)
#----------------------------------------------------------------------------------------------------------------------#

def dumpOutput(input, output, fileName):
    dumpObj = {}
    dumpObj['input'] = input
    dumpObj['output'] = output
    with open(fileName, 'wb') as outfile:
        pickle.dump(dumpObj, outfile)