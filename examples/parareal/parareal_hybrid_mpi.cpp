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
#include <pfasst/encap/vector.hpp>

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
          
          shared_ptr<Encapsulation<double>> u; // vector with current numerical solution at fine level
          vector<shared_ptr<Encapsulation<double>>> uExact; // vector with exact solution at fine level
          shared_ptr<Encapsulation<double>> err;
          error_map errors;
          
          double dt; // time step
          size_t nfineiters; // number of sdc iterations for fine propagator
          size_t ncrseiters; // number of sdc iterations for coarse propagator
          
          size_t numTiters; // number of parareal repetitions
          
          double timeMeasure;
          
          shared_ptr<SDC<time>> getSDC(const size_t nnodes, 
                                       const pfasst::quadrature::QuadratureType quad_type, 
                                       const size_t ndofs, const double abs_res_tol, 
                                       const double rel_res_tol, bool coarse)
          {
            auto sdc = make_shared<SDC<time>>();
            
            auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
            auto factory = make_shared<VectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);
            sweeper->set_quadrature(quad);
            sweeper->set_factory(factory);
            sweeper->set_controller(&*sdc);
            sweeper->setup(coarse);
            
            sdc->add_level(sweeper);
            
            if(commRank == 0) {
              auto q0 = sweeper->get_start_state();
              sweeper->exact(q0, 0.0);
            }
            
            return sdc;
          }
          
          void init(const double t_end, const double dt, 
                    const size_t nnodes, const size_t nnodesCoarse,
                    const pfasst::quadrature::QuadratureType quad_type,
                    const size_t ndofs_fine, const size_t ndofs_coarse,
                    const size_t nfineiters, const size_t ncrseiters,
                    const double abs_res_tol, const double rel_res_tol)
          {
            this->dt = dt;
            this->nfineiters = nfineiters;
            this->ncrseiters = ncrseiters;
            
            size_t numTSlices = t_end/dt;
            numTiters = numTSlices/commSize;
            if(numTSlices % commSize != 0) numTiters++;
            
            if(commRank == 0) {
              CLOG(INFO, "Parareal") << "numTSlices: " << numTSlices;
              CLOG(INFO, "Parareal") << "numTiters: " << numTiters;
            }
            
            factory = make_shared<VectorFactory<double>>(ndofs_fine);
            
            fine = getSDC(nnodes, quad_type, ndofs_fine, abs_res_tol, rel_res_tol, false);
            coarse = getSDC(nnodesCoarse, quad_type, ndofs_coarse, abs_res_tol, rel_res_tol, true);
            
            err = factory->create(solution);
            transfer = make_shared<SpectralTransfer1D<>>();
            
            u = factory->create(solution);
            
            uExact.clear();
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs_fine);
            for(size_t i = 0; i < numTiters; i++) {
              uExact.push_back(factory->create(solution));
              sweeper->exact(uExact[i], (i*commSize + commRank + 1)*dt);
            }
          }
          
          void do_coarse(shared_ptr<Encapsulation<>> start_state,
                         shared_ptr<Encapsulation<>> end_state, bool predict) 
          {
            coarse->set_step(0);
            coarse->set_iteration(predict ? 0 : 1);
            EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
            
            if(start_state) {
              transfer->restrict(coarseSweeper.get_start_state(), start_state);
              coarseSweeper.reevaluate();
            }
            
            if(predict) CLOG(INFO, "Parareal") << "Predict";
            else CLOG(INFO, "Parareal") << "Coarse-Sweep";
            
            coarse->run();
            
            if(as_vector<double>(coarseSweeper.get_end_state()).size() == as_vector<double>(end_state).size()) {
              end_state->copy(coarseSweeper.get_end_state());
            } 
            else {
              transfer->interpolate(end_state, coarseSweeper.get_end_state());
            }
          }
          
          void do_fine(shared_ptr<Encapsulation<>> start_state,
                       shared_ptr<Encapsulation<>> end_state) 
          {
            fine->set_step(0);
            fine->set_iteration(1);
            EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
            
            if(start_state) {
              fineSweeper.get_start_state()->copy(start_state);
              fineSweeper.reevaluate();
            }
                
            CLOG(INFO, "Parareal") << "Fine-Sweep";
            
            fine->run();
            
            end_state->copy(fineSweeper.get_end_state());
          }
          
          shared_ptr<Encapsulation<double>> recvState(const size_t k, const size_t j, bool* prec_done)
          {
            int t = tag(k, j, commRank);
            CLOG(INFO, "Parareal") << "Recv tag: " << t;
            shared_ptr<Encapsulation<double>> state = factory->create(solution);
            state->recv(comm, t, true);
            
            if(k > 0 && commRank > 0)
            {
              comm->status->recv(t);
              *prec_done = comm->status->get_converged(commRank - 1);
            }
            
            return state;
          }
          
          int tag(const size_t k, const size_t j, int commRank) 
          {
            return k * 10000 + commSize * j + commRank;
          }
          
          void echo_error(const size_t k, const size_t j)
          {
            size_t n_global = commSize * j + commRank;
            err->copy(uExact[j]);
            err->saxpy(-1.0, u);
            errors.insert(vtype(ktype(n_global, k-1), err->norm0()));
            CLOG(INFO, "Parareal") << "Error: " << err->norm0() << " k: " << k << " n: " << n_global;
          }
          
        public:
          
          static void init_opts()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts();
            pfasst::config::options::add_option<size_t>("Parareal", "num_par_iter", "Number of Parareal iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "spatial_dofs_coarse", "Number of spatial degrees of freedom at coarse level");
            pfasst::config::options::add_option<size_t>("Parareal", "num_nodes_coarse", "Number of collocation nodes for coarse sweeper");
            pfasst::config::options::add_option<size_t>("Parareal", "num_fine_iter", "Number of Fine Sweep iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "num_crse_iter", "Number of Coarse Sweep iterations");
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
          
          void run_parareal(const size_t npariters, const size_t nfineiters, 
                            const size_t ncrseiters, const double dt, const double t_end,
                            const size_t nnodes, const size_t nnodesCoarse,
                            const pfasst::quadrature::QuadratureType quad_type,
                            const size_t ndofs_fine, const size_t ndofs_coarse,
                            const double abs_res_tol, const double rel_res_tol)
          {
            if(!comm || commSize == 1) {
              CLOG(ERROR, "Parareal") << "Number of processors must be greater than 1";
              MPI_Finalize();
              exit(-1);
            }
            
            init(t_end, dt, nnodes, nnodesCoarse,
                 quad_type, ndofs_fine, ndofs_coarse, nfineiters, 
                 ncrseiters, abs_res_tol, rel_res_tol);
            
            auto coarseState = factory->create(solution);
            auto fineState = factory->create(solution);
            auto diff = factory->create(solution);
            shared_ptr<Encapsulation<double>> startState;
            
            double res;
            
            MPI_Barrier(MPI_COMM_WORLD);
            timeMeasure = MPI_Wtime();
            
            bool prec_done = false;
            bool done = false;
            
            for(size_t j = 0; j < numTiters; j++) { // loop over time blocks
              
              size_t nglobal = commSize * j + commRank;
              
              // set the time step
              coarse->set_duration(nglobal*dt, nglobal*dt+dt, dt, 2);
              fine->set_duration(nglobal*dt, nglobal*dt+dt, dt, 2);
              
              for(size_t k = 0; k < npariters && !done; k++) { // loop over parareal iterations
                
                if(k > 0) {
                  if(k == 1) transfer->PolyInterpMixin<time>::interpolate(fine->get_finest(), coarse->get_finest(), true);
                  do_fine(startState, fineState);
                  transfer->PolyInterpMixin<time>::restrict(coarse->get_finest(), fine->get_finest(), true);
                  diff->copy(u);
                  
//                   auto& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
//                   auto& crseSweeper = as_encap_sweeper<time>(coarse->get_finest());
//                   size_t nfine = fineSweeper.get_nodes().size();
//                   size_t ncrse = crseSweeper.get_nodes().size();
//                   CLOG(INFO, "Parareal") << "after fine sweep";
//                   for (size_t m = 0; m < nfine; m++) CLOG(INFO, "Parareal") << ((vector<double>)(as_vector<double>(fineSweeper.get_state(m))));
//                   CLOG(INFO, "Parareal") << "restricted";
//                   for (size_t m = 0; m < ncrse; m++) CLOG(INFO, "Parareal") << ((vector<double>)(as_vector<double>(crseSweeper.get_state(m))));
                }
                
                if(!prec_done && (commRank > 0 || (k == 0 && j > 0))) startState = recvState(k,j, &prec_done);
                
                coarse->get_finest()->save();
                if(k == 0) coarse->set_duration(nglobal*dt, nglobal*dt+dt, dt, 1);
                do_coarse(startState, coarseState, k == 0);
                if(k == 0) coarse->set_duration(nglobal*dt, nglobal*dt+dt, dt, 2);
                u->copy(coarseState);
                
                echo_error(k, j);
                
                if(k > 0) {
                  // calc the residium
                  diff->saxpy(-1.0, u);
                  res = diff->norm0();
                  done = res < abs_res_tol;
                  
                  CLOG(INFO, "Parareal") << "Residium: " << res;
                  CLOG(INFO, "Parareal") << "Done: " << done;
                }
                
                if(commRank < commSize - 1 && (nglobal+1)*dt <= t_end) {
                  int t = tag(k, j, commRank+1);
                  CLOG(INFO, "Parareal") << "Send tag: " << t;
                  u->send(comm, t, true);
                  
                  if(k > 0) {
                    comm->status->set_converged(done);
                    comm->status->send(t);
                  }
                }
              } // loop over parareal iterations
              
              comm->status->clear();
              prec_done = false;
              done = false; 
              
              if(j < numTiters - 1 && commRank == commSize - 1) {
                CLOG(INFO, "Parareal") << "Send tag: " << tag(0, j+1, 0);
                u->send(comm, tag(0, j+1, 0), true); // Send to proc with rank 0
              }
            } // loop over time blocks
              
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
  
  pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.set_comm(&comm);
    
  const size_t  npariters = pfasst::config::get_value<size_t>("num_par_iter", 20);
  const size_t  nfineiters = pfasst::config::get_value<size_t>("num_fine_iter", 10);
  const size_t  ncrseiters = pfasst::config::get_value<size_t>("num_crse_iter", 10);
  const size_t  nnodes = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t  nnodesCoarse = pfasst::config::get_value<size_t>("num_nodes_coarse", 5);
  const size_t  ndofs_fine  = pfasst::config::get_value<size_t>("spatial_dofs", 64);
  const size_t  ndofs_coarse  = pfasst::config::get_value<size_t>("spatial_dofs_coarse", 64);
  auto const quad_type = \
    pfasst::config::get_value<pfasst::quadrature::QuadratureType>("nodes_type", pfasst::quadrature::QuadratureType::GaussLobatto);
  
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);
  
  parareal.run_parareal(npariters, nfineiters, ncrseiters, dt, 
                        t_end, nnodes, nnodesCoarse,
                        quad_type, ndofs_fine, ndofs_coarse,
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
