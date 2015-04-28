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
          
          double dt; // time step
          size_t nfineiters; // number of sdc iterations for fine propagator
          size_t ncrseiters; // number of sdc iterations for coarse propagator
          
          size_t numTiters; // number of parareal repetitions
          size_t ndt_per_proc;
          
          double timeMeasure;
          
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
          
          void init(const size_t ndt_per_proc, const double t_end, 
                    const double dt, const size_t nnodes,
                    const pfasst::quadrature::QuadratureType quad_type,
                    const size_t ndofs_fine, const size_t ndofs_coarse,
                    const size_t nfineiters, const size_t ncrseiters,
                    const double abs_res_tol, const double rel_res_tol)
          {
            this->ndt_per_proc = ndt_per_proc;
            this->dt = dt;
            this->nfineiters = nfineiters;
            this->ncrseiters = ncrseiters;
            
            size_t numTSlices = t_end/dt;
            size_t numBlocks = numTSlices/ndt_per_proc;
            if(numTSlices % ndt_per_proc != 0) numBlocks++;
            numTiters = numBlocks/commSize;
            if(numBlocks % commSize != 0) numTiters++;
            
            if(commRank == 0) {
              CLOG(INFO, "Parareal") << "numTSlices: " << numTSlices;
              CLOG(INFO, "Parareal") << "numBlocks: " << numBlocks;
              CLOG(INFO, "Parareal") << "numTiters: " << numTiters;
            }
            
            factory = make_shared<MPIVectorFactory<double>>(ndofs_fine);
            
            fine = getSDC(nnodes, quad_type, ndofs_fine, abs_res_tol, rel_res_tol);
            coarse = getSDC((nnodes + 1) / 2, quad_type, ndofs_coarse, abs_res_tol, rel_res_tol);
            
            err = factory->create(solution);
            transfer = make_shared<SpectralTransfer1D<>>();
            
            u.clear();
            uFine.clear();
            uCoarse.clear();
            uExact.clear();
            
            for(size_t i = 0 ; i < ndt_per_proc; i++) {
              u.push_back(factory->create(solution));
              uFine.push_back(factory->create(solution));
              uCoarse.push_back(factory->create(solution));
            }
            
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs_fine);
            for(size_t i = 0; i < numTiters; i++) {
              for(size_t j = 0; j < ndt_per_proc; j++) {
                uExact.push_back(factory->create(solution));
                sweeper->exact(uExact[i*ndt_per_proc+j], ((i*commSize+commRank)*ndt_per_proc + j + 1)*dt);
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
          
          void do_fine(shared_ptr<Encapsulation<>> start_state,
                       shared_ptr<Encapsulation<>> end_state,
                       const double startTime) 
          {
            fine->set_duration(startTime, startTime+dt, dt, nfineiters);
            EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
            
            if(start_state) fineSweeper.get_start_state()->copy(start_state);
                
            CLOG(INFO, "Parareal") << "Fine-Sweep";
            
            fine->run();
            
            end_state->copy(fineSweeper.get_end_state());
          }
          
          shared_ptr<Encapsulation<double>> recvState(const size_t k, const size_t j, bool* prec_done)
          {
            CLOG(INFO, "Parareal") << "Recv tag: " << tag(k, j, commRank);
            shared_ptr<Encapsulation<double>> state = factory->create(solution);
            CLOG(INFO, "Parareal") << "passedTime-Recv mpi_vector before: " << MPI_Wtime() - timeMeasure;
            state->recv(comm, tag(k, j, commRank), true);
            CLOG(INFO, "Parareal") << "passedTime-Recv mpi_vector after: " << MPI_Wtime() - timeMeasure;
            
            if(commRank > 0)
            {
              comm->status->recv();
              *prec_done = comm->status->get_converged(commRank - 1);
            }
            
            return state;
          }
          
          int tag(const size_t k, const size_t j, int commRank) 
          {
            return k * 1000 + (commSize * j + commRank) * ndt_per_proc + 10;
          }
          
          size_t get_n_local(const size_t n) 
          {
            return commRank * ndt_per_proc + n;
          }
          
          size_t get_n_global(const size_t n, const size_t j) 
          {
            return commSize * ndt_per_proc * j + get_n_local(n);
          }
          
          void echo_error(const size_t n, const size_t k, const size_t j) 
          {
            size_t n_global = get_n_global(n, j);
            err->copy(uExact[j * ndt_per_proc + n]);
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
          
          void run_parareal(const size_t ndt_per_proc, const size_t npariters,
                            const size_t nfineiters, const size_t ncrseiters,
                            const double dt, const double t_end, const size_t nnodes,
                            const pfasst::quadrature::QuadratureType quad_type,
                            const size_t ndofs_fine, const size_t ndofs_coarse,
                            const double abs_res_tol, const double rel_res_tol)
          {
            init(ndt_per_proc, t_end, dt, nnodes, quad_type, 
                 ndofs_fine, ndofs_coarse, nfineiters, 
                 ncrseiters, abs_res_tol, rel_res_tol);
            
            if(!comm || commSize == 1) {
              CLOG(ERROR, "Parareal") << "Number of processors must be greater than 1";
              MPI_Finalize();
              exit(-1);
            }
            
            shared_ptr<Encapsulation<double>> coarseState = factory->create(solution);
            shared_ptr<Encapsulation<double>> diff = factory->create(solution);
            shared_ptr<Encapsulation<double>> startState;
            shared_ptr<Encapsulation<double>> delta = factory->create(solution);
            
            double res;
            
            MPI_Barrier(MPI_COMM_WORLD);
            timeMeasure = MPI_Wtime();
            
            bool prec_done = false;
            bool done = false;
            
            size_t nstart = 0;
            
            for(size_t j = 0; j < numTiters; j++) { // loop over time blocks
              
              for(size_t k = 0; k < npariters && !done; k++) { // loop over parareal iterations
                
                if(prec_done) nstart++;
                if(nstart > ndt_per_proc || get_n_global(0, j)*dt >= t_end) break;
                
                if(k > 0) {
                  for(size_t n = nstart; n < ndt_per_proc; n++) {
                    
                    if(get_n_local(n - nstart) < k-1) continue;
                    
                    do_fine(n > 0 ? u[n-1] : startState, uFine[n], get_n_global(n, j)*dt);
                  }
                }
                
                if(get_n_local(0) >= k && ((commRank > 0 && !prec_done) || (j > 0 && k == 0))) {
                  startState = recvState(k,j, &prec_done);
                }
                
                if(k > 0) {
                  diff->copy(u[ndt_per_proc - 1]); // calculate residium for break condition
                }
                  
                for(size_t n = nstart; n < ndt_per_proc; n++) {
                  
                  if(k > 0 && get_n_local(n - nstart) < k-1) continue;
                  
                  if(get_n_local(n - nstart) >= k) {
                    
                    do_coarse(n > 0 ? u[n-1] : startState, coarseState, get_n_global(n, j)*dt);
                    
                    if(k == 0) {
                      u[n]->copy(coarseState);
                    }
                    else {
                      u[n]->copy(uFine[n]);
                      
                      delta->copy(coarseState);
                      delta->saxpy(-1.0, uCoarse[n]);
                      
                      // apply delta correction
                      u[n]->saxpy(1.0, delta);
                      
                      CLOG(INFO, "Parareal") << "Delta-Norm: " << delta->norm0();
                    }
                    uCoarse[n]->copy(coarseState);  
                  }
                  else {
                    u[n]->copy(uFine[n]);
                  }
                  echo_error(n, k, j);
                }
                
                if(k > 0) {
                  // calc the residium
                  diff->saxpy(-1.0, u[ndt_per_proc - 1]);
                  res = diff->norm0();
                  done = res < abs_res_tol;
                  
                  CLOG(INFO, "Parareal") << "Residium: " << res;
                  CLOG(INFO, "Parareal") << "Done: " << done;
                }
                
                if(commRank < commSize - 1 && (get_n_global(ndt_per_proc - 1, j)+1)*dt < t_end) {
                  CLOG(INFO, "Parareal") << "Send tag: " << tag(k, j, commRank+1);
                  u[ndt_per_proc-1]->send(comm, tag(k, j, commRank+1), false);
                  
                  done = done || get_n_local(ndt_per_proc - 1 - nstart) < k;
                  comm->status->set_converged(done);
                  comm->status->send();
                  CLOG(INFO, "Parareal") << "passedTime-Send: " << MPI_Wtime() - timeMeasure;
                }
              } // loop over parareal iterations
              
              comm->status->clear();
              prec_done = false;
              done = false;
              nstart = 0;
              
              if(j < numTiters - 1 && commRank == commSize - 1) {
                CLOG(INFO, "Parareal") << "Send tag: " << tag(0, j+1, 0);
                u[ndt_per_proc-1]->send(comm, tag(0, j+1, 0), false); // Send to proc with rank 0
                CLOG(INFO, "Parareal") << "passedTime-Send: " << MPI_Wtime() - timeMeasure;
              }
            } // loop over ring-blocking iterations
            
            timeMeasure = MPI_Wtime() - timeMeasure;
            CLOG(INFO, "Parareal") << "time Measurement: " << timeMeasure;
          } // function run_parareal
      }; // Class Parareal
    } // namespace parareal
  } // namespace examples
} // namespace pfasst


void run_example()
{
  const double  dt     = pfasst::config::get_value<double>("dt", 0.01);
  const double  t_end     = pfasst::config::get_value<double>("tend", 1.0);
  
  // dt must be a divisor of t_end
  if(t_end - (int)(t_end/dt)*dt != 0) {
    MPI_Finalize();
    CLOG(INFO, "Parareal") << "Input-error: dt must be a divisor of t_end";
    exit(-1);
  }
    
  const size_t  npariters = pfasst::config::get_value<size_t>("num_par_iter", 5);
  const size_t  nfineiters = pfasst::config::get_value<size_t>("num_fine_iter", 10);
  const size_t  ncrseiters = pfasst::config::get_value<size_t>("num_crse_iter", 10);
  const size_t  nnodes = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t  ndofs_fine  = pfasst::config::get_value<size_t>("spatial_dofs", 64);
  const size_t  ndofs_coarse  = pfasst::config::get_value<size_t>("spatial_dofs_coarse", 64);
  const size_t  ndt_per_proc  = pfasst::config::get_value<size_t>("ndt_per_proc", 1);
  auto const quad_type = \
    pfasst::config::get_value<pfasst::quadrature::QuadratureType>("nodes_type", pfasst::quadrature::QuadratureType::GaussLobatto);
  
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 1e-10);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);
    
  pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.set_comm(&comm);
  parareal.run_parareal(ndt_per_proc, npariters, nfineiters, ncrseiters, 
                        dt, t_end, nnodes, quad_type, ndofs_fine, ndofs_coarse,
                        abs_res_tol, rel_res_tol);
}


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  
  pfasst::init(argc, argv,
               pfasst::examples::parareal::Parareal<>::init_opts,
               pfasst::examples::parareal::Parareal<>::init_logs);
  
  run_example();
  
  fftw_cleanup();
  MPI_Finalize();
}
#endif
