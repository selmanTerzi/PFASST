#include <iostream>
#include <memory>
#include <vector>
using namespace std;

#include <mpi.h>
#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/encap_sweeper.hpp>
#include <pfasst/mpi_communicator.hpp>
#include <pfasst/encap/mpi_vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

using namespace pfasst::encap;

namespace pfasst 
{
  namespace examples 
  {
    namespace parareal 
    {
      template<typename time = time_precision>
      class Parareal
      {
        
        private:
          
          ICommunicator* comm;
          
          int commSize;
          int commRank;
          
          shared_ptr<SDC<time>> coarse; // coarse-level controller
          shared_ptr<SDC<time>> fine; // fine-level controller
          
          shared_ptr<SpectralTransfer1D<time>> transfer;
          
          shared_ptr<EncapFactory<>> factory; // fine level factory
          
          vector<shared_ptr<Encapsulation<double>>> u; // vector with current numerical solution at fine level
          vector<shared_ptr<Encapsulation<double>>> uCoarse; // vector with interpolated coarse-Sweep results
          vector<shared_ptr<Encapsulation<double>>> uFine; // vector with fine-Sweep results
          vector<shared_ptr<Encapsulation<double>>> uExact; // vector with exact solution at fine level
          shared_ptr<Encapsulation<double>> err;
          error_map errors;
          
          shared_ptr<Encapsulation<double>> startState; // start_state for current time-slice
          vector<bool> converged;
          vector<bool> done;
          
          double dt; // time step
          size_t nfineiters; // number of sdc iterations for fine propagator
          size_t ncrseiters; // number of sdc iterations for coarse propagator
          
          size_t n_global; // global number of the current time slice
          size_t numTiters; // number of parareal repetitions
          
          shared_ptr<SDC<time>> getSDC(const size_t nnodes, 
                                       const pfasst::quadrature::QuadratureType quad_type, 
                                       const size_t ndofs, const double abs_res_tol, 
                                       const double rel_res_tol)
          {
            auto sdc = make_shared<SDC<time>>();
            
            auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
            auto factory = make_shared<MPIVectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);
            sweeper->set_quadrature(quad);
            sweeper->set_factory(factory);
            
            sdc->add_level(sweeper);
            sdc->setup();
            
            if(commRank == 0) {
              auto q0 = sweeper->get_start_state();
              sweeper->exact(q0, 0.0);
            }
            
            return sdc;
          }
          
          void init(const size_t nsteps, const double t_end, 
                    const double dt, const size_t nnodes,
                    const pfasst::quadrature::QuadratureType quad_type,
                    const size_t ndofs_fine, const size_t ndofs_coarse,
                    const size_t nfineiters, const size_t ncrseiters,
                    const double abs_res_tol, const double rel_res_tol)
          {
            
            size_t numTSlices = t_end/dt;
            size_t numBlocks = numTSlices/nsteps;
            if(numTSlices % nsteps != 0) numBlocks++;
            numTiters = numBlocks/commSize;
            if(numBlocks % commSize != 0) numTiters++;
            
            if(commRank == 0) {
              CLOG(INFO, "Parareal") << "numTSlices: " << numTSlices;
              CLOG(INFO, "Parareal") << "numBlocks: " << numBlocks;
              CLOG(INFO, "Parareal") << "numTiters: " << numTiters;
            }
            
            factory = make_shared<MPIVectorFactory<double>>(ndofs_fine);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs_fine);
            this->dt = dt;
            this->nfineiters = nfineiters;
            this->ncrseiters = ncrseiters;
            
            fine = getSDC(nnodes, quad_type, ndofs_fine, abs_res_tol, rel_res_tol);
            coarse = getSDC(nnodes, quad_type, ndofs_coarse, abs_res_tol, rel_res_tol);
            
            err = factory->create(solution);
            transfer = make_shared<SpectralTransfer1D<>>();
            
            u.clear();
            uFine.clear();
            uCoarse.clear();
            uExact.clear();
            
            for(size_t i = 0 ; i < nsteps; i++) {
              u.push_back(factory->create(solution));
              uFine.push_back(factory->create(solution));
              uCoarse.push_back(factory->create(solution));
              converged.push_back(false);
              done.push_back(false);
            }
            
            for(size_t i = 0; i < numTiters; i++) {
              for(size_t j = 0; j < nsteps; j++) {
                cout << i << " " << j << endl;
                uExact.push_back(factory->create(solution));
                sweeper->exact(uExact[i*nsteps+j], ((i*commSize+commRank)*nsteps + j + 1)*dt);
              }
            }
          }
          
          void do_coarse(shared_ptr<Encapsulation<>> start_state,
                         shared_ptr<Encapsulation<>> end_state,
                         const double startTime) 
          {
            coarse->set_duration(startTime, startTime+dt, dt, ncrseiters);
            EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
            
            if(start_state) transfer->restrict(coarseSweeper.get_start_state(), start_state);
                
            CLOG(INFO, "Parareal") << "Coarse-Sweep";
            
            coarse->run();
            
            if(as_vector<double>(coarseSweeper.get_end_state()).size() == as_vector<double>(end_state).size()) {
              end_state->copy(coarseSweeper.get_end_state());
            } 
            else {
              transfer->interpolate(end_state, coarseSweeper.get_end_state());
            }
          }
          
          bool do_fine(shared_ptr<Encapsulation<>> start_state,
                       shared_ptr<Encapsulation<>> end_state,
                       const double startTime) 
          {
            fine->set_duration(startTime, startTime+dt, dt, nfineiters);
            EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
            
            if(start_state) fineSweeper.get_start_state()->copy(start_state);
                
            CLOG(INFO, "Parareal") << "Fine-Sweep";
            fine->run();
            
            end_state->copy(fineSweeper.get_end_state());
            return fineSweeper.converged();
          }
          
          void get_state(const size_t n, const size_t k, const size_t j)
          {
            if((commRank > 0 || n > 0 || j > 0) && !startState) startState = factory->create(solution); 
            
            if(n > 0) {
              startState->copy(u[n-1]);
            }
            else if(j > 0 || commRank > 0) {
              startState->recv(comm, k * 1000 + n_global, true);
            }
          }
          
          void echo_error(size_t n, size_t nsteps, size_t k, size_t j) {
            err->copy(uExact[j * nsteps + n]);
            err->saxpy(-1.0, u[n]);
            errors.insert(vtype(ktype(n_global, k-1), err->norm0()));
            CLOG(INFO, "Parareal") << "Error: " << err->norm0() << " k: " << k << " n: " << n_global;
          }
          
        public:
          
          static void init_opts()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts();
            pfasst::config::options::add_option<size_t>("Parareal", "num_par_iter", "Number of Parareal iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "spatial_dofs_coarse", "Number of spatial degrees of freedom at coarse level");
            pfasst::config::options::add_option<size_t>("Parareal", "num_fine_iter", "Number of Fine Sweep iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "num_crse_iter", "Number of Coarse Sweep iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "ndt_per_proc", "Number of time slices per process");
            pfasst::config::options::add_option<time>("Parareal", "t_end", "End time");
          }
          
          static void init_logs()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_logs();
            pfasst::log::add_custom_logger("Parareal");
          }
          
          void set_comm(ICommunicator* comm)
          {
            this->comm = comm;
            commRank = comm->rank();
            commSize = comm->size();
          }
          
          void run_parareal(const size_t nsteps, const size_t niters,
                            const size_t nfineiters, const size_t ncrseiters,
                            const double dt, const double t_end, const size_t nnodes,
                            const pfasst::quadrature::QuadratureType quad_type,
                            const size_t ndofs_fine, const size_t ndofs_coarse,
                            const double abs_res_tol, const double rel_res_tol)
          {
            init(nsteps, t_end, dt, nnodes, quad_type, 
                 ndofs_fine, ndofs_coarse, nfineiters, 
                 ncrseiters, abs_res_tol, rel_res_tol);
            
            if(!comm || commSize == 1) {
              cout << "Number of processors must be greater than 1";
              MPI_Finalize();
              exit(-1);
            }
            
            shared_ptr<Encapsulation<double>> coarseState = factory->create(solution);
            shared_ptr<Encapsulation<double>> fineState = factory->create(solution);
            shared_ptr<Encapsulation<double>> delta = factory->create(solution);
            
            double timeMeasure = MPI_Wtime();
            
            for(size_t j = 0; j < numTiters; j++) {
              
              for(size_t k = 0; k < niters; k++) {
                
                if(done[nsteps - 1]) break;
                
                for(size_t n = 0; n < nsteps; n++) {
                  
                  size_t n_local = commRank * nsteps + n;
                  n_global = j * commSize * nsteps + n_local;
                  
                  if(n_global * dt > t_end) break; // break if t_end is reached
                  if(done[n] || k > n_global + 1) continue;
                  
                  CLOG(INFO, "Parareal") << " k: " << k << " n: " << n_global << " j: " << j;
                  
                  if(n_local >= k) {
                    if(commRank == 0 || comm->status->previous_is_iterating()) {
                      get_state(n, k, j);
                    }
                    do_coarse(startState, coarseState, n_global*dt);
                    if(!converged[n] && k < niters - 1) converged[n] = do_fine(startState, fineState, n_global*dt);
                  }
                    
                  if(k > 0) {
                    u[n]->copy(uFine[n]);
                    
                    if(n >= k) {
                      
                      delta->copy(coarseState);
                      delta->saxpy(-1.0, uCoarse[n]);
                      
                      CLOG(INFO, "Parareal") << "delta Norm: " << delta->norm0() << " k: " << k << " n: " << n_global;
                      
                      // apply delta correction
                      u[n]->saxpy(1.0, delta);
                    }
                    
                    if(converged[n]) {
                      done[n] = true;
                      if(n == nsteps-1) comm->status->set_converged(true);
                    }
                  }
                  else {
                    u[n]->copy(coarseState);
                  }
                  
                  if(n_local >= k) {
                    uCoarse[n]->copy(coarseState);
                    uFine[n]->copy(fineState);
                  }
                  
                  echo_error(n, nsteps, k, j);
                } // loop over time slices
                
                if(commRank < commSize - 1) {
                  u[nsteps-1]->send(comm, k * 1000 + n_global+1, false);
                }
              } // loop over parareal iterations
              
              comm->status->clear();
              
              if(j < numTiters - 1) {
                for(size_t i = 0; i < nsteps; i++) {
                  converged[i] = done[i] = false;
                }
                if(commRank == commSize - 1) {
                  u[nsteps-1]->send(comm, n_global+1, false); // Send to proc with rank 0
                }    
              }
            } // loop over ring-blocking iterations
            
            timeMeasure = MPI_Wtime() - timeMeasure;
            CLOG(INFO, "Parareal") << "time Measurement: " << timeMeasure;
          } // function run_parareal
      }; // Class Parareal
    } // namespace parareal
  } // namespace examples
} // namespace pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  
  pfasst::init(argc, argv,
               pfasst::examples::parareal::Parareal<>::init_opts,
               pfasst::examples::parareal::Parareal<>::init_logs);
  
  const double  dt     = pfasst::config::get_value<double>("dt", 0.01);
  const double  t_end     = pfasst::config::get_value<double>("t_end", 1.0);
  
  // dt must be a divisor of t_end
  if(t_end - (int)(t_end/dt)*dt != 0) {
    MPI_Finalize();
    cout << t_end - (int)(t_end/dt)*dt << endl;
    cout << "dt must be a divisor of t_end" << endl;
    exit(-1);
  }
    
  const size_t  npariters = pfasst::config::get_value<size_t>("num_par_iter", 5);
  const size_t  nfineiters = pfasst::config::get_value<size_t>("num_fine_iter", 10);
  const size_t  ncrseiters = pfasst::config::get_value<size_t>("num_crse_iter", 10);
  const size_t  nnodes = pfasst::config::get_value<size_t>("num_nodes", 3);
  const size_t  ndofs_fine  = pfasst::config::get_value<size_t>("spatial_dofs", 64);
  const size_t  ndofs_coarse  = pfasst::config::get_value<size_t>("spatial_dofs_coarse", 64);
  const size_t  ndt_per_proc  = pfasst::config::get_value<size_t>("ndt_per_proc", 1);
  auto const quad_type = \
    pfasst::config::get_value<pfasst::quadrature::QuadratureType>("nodes_type", pfasst::quadrature::QuadratureType::GaussLobatto);
  
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);
    
  pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.set_comm(&comm);
  parareal.run_parareal(ndt_per_proc, npariters, nfineiters, ncrseiters, 
                        dt, t_end, nnodes, quad_type, ndofs_fine, ndofs_coarse,
                        abs_res_tol, rel_res_tol);
  
  fftw_cleanup();
  MPI_Finalize();
}
#endif
