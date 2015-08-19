#include <memory>
#include <vector>
#include <ctime>
using namespace std;

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/encap_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"
#include <fftw3.h>

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
          shared_ptr<SDC<time>> coarse; // coarse-level controller
          shared_ptr<SDC<time>> fine; // fine-level controller
          
          shared_ptr<SpectralTransfer1D<time>> transfer;
          
          shared_ptr<EncapFactory<>> factory; // fine level factory
          
          vector<shared_ptr<Encapsulation<double>>> u; // vector with current numerical solution at fine level
          vector<shared_ptr<Encapsulation<double>>> uCoarse; // vector with interpolated coarse-Sweep results
          vector<shared_ptr<Encapsulation<double>>> uFine; // vector with fine-Sweep results
          vector<shared_ptr<Encapsulation<double>>> uExact; // vector with exact solution at fine level
          
          shared_ptr<Encapsulation<double>> err;
          
          vector<bool> done;
          
          double dt; // time step
          size_t nfineiters; // number of sdc iterations for fine propagator
          size_t ncrseiters; // number of sdc iterations for coarse propagator
          
          shared_ptr<SDC<time>> getSDC(const size_t nnodes,
                                       const quadrature::QuadratureType quad_type,
                                       const size_t ndofs, 
                                       const double abs_res_tol, const double rel_res_tol)
          {
            auto sdc = make_shared<SDC<time>>();
            auto quad = quadrature::quadrature_factory(nnodes, quad_type);
            auto factory = make_shared<VectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);
            sweeper->set_quadrature(quad);
            sweeper->set_factory(factory);
            
            sdc->add_level(sweeper);
            sdc->setup();
            
            auto q0 = sweeper->get_start_state();
            sweeper->exact(q0, 0.0);
            
            return sdc;
          }
        
          void init(const double tend, const double dt,
                    const size_t nnodes, const size_t nnodesCoarse,
                    const quadrature::QuadratureType quad_type,
                    const size_t ndofs_fine, const size_t ndofs_coarse,
                    const size_t nfineiters, const size_t ncrseiters,
                    const double abs_res_tol, const double rel_res_tol)
          {
            
            factory = make_shared<VectorFactory<double>>(ndofs_fine);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs_fine);
            this->dt = dt;
            this->nfineiters = nfineiters;
            this->ncrseiters = ncrseiters;
            
            fine = getSDC(nnodes, quad_type, ndofs_fine, abs_res_tol, rel_res_tol);
            coarse = getSDC(nnodesCoarse, quad_type, ndofs_coarse, abs_res_tol, rel_res_tol);
            
            err = factory->create(solution);
            transfer = make_shared<SpectralTransfer1D<>>();
            
            u.clear();
            uFine.clear();
            uCoarse.clear();
            uExact.clear();
            
            size_t nsteps = tend/dt;
            
            for(size_t i = 0 ; i < nsteps; i++) {
              u.push_back(factory->create(solution));
              uFine.push_back(factory->create(solution));
              uCoarse.push_back(factory->create(solution));
              uExact.push_back(factory->create(solution));
              sweeper->exact(uExact[i],(i+1)*dt);
            }
          }
          
          void do_coarse(shared_ptr<Encapsulation<>> start_state,
                         shared_ptr<Encapsulation<>> end_state,
                         const double startTime) 
          {
            coarse->set_duration(startTime, startTime+dt, dt, ncrseiters);
            EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
            
            if(start_state) transfer->restrict(coarseSweeper.get_start_state(), start_state);
                
            CLOG(INFO, "Parareal")  << "Coarse-Sweep";
            coarse->run();
            
            transfer->interpolate(end_state, coarseSweeper.get_end_state());
          }
          
          void do_fine(shared_ptr<Encapsulation<>> start_state,
                       shared_ptr<Encapsulation<>> end_state,
                       const double startTime) 
          {
            fine->set_duration(startTime, startTime+dt, dt, nfineiters);
            EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
            
            if(start_state) fineSweeper.get_start_state()->copy(start_state);
                
            CLOG(INFO, "Parareal")  << "Fine-Sweep";
            fine->run();
            
            end_state->copy(fineSweeper.get_end_state());
          }
          
          void echo_error(const size_t n, const size_t k, const double residual)
          {
            err->copy(uExact[n]);
            err->saxpy(-1.0, u[n]);
            CLOG(INFO, "Parareal") << "ErrorParareal: step: " << n+1<< " iter: " << k << " residual: " << residual << " err: " << err->norm0();
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
        
          void run_parareal(const double tend, const size_t niters,
                            const size_t nfineiters, const size_t ncrseiters,
                            const double dt, const size_t nnodes, const size_t nnodesCoarse,
                            const quadrature::QuadratureType quad_type,
                            const size_t ndofs_fine, const size_t ndofs_coarse,
                            const double abs_res_tol, const double rel_res_tol)
          {
            
            init(tend, dt, nnodes, nnodesCoarse, quad_type, ndofs_fine, ndofs_coarse, 
                 nfineiters, ncrseiters, abs_res_tol, rel_res_tol);
            
            shared_ptr<Encapsulation<double>> coarseState = factory->create(solution);
            shared_ptr<Encapsulation<double>> fineState = factory->create(solution);
            shared_ptr<Encapsulation<double>> delta = factory->create(solution);
            shared_ptr<Encapsulation<double>> diff = factory->create(solution);
            double res = 0.0;
            
            clock_t timeMeasure = clock();
            
            size_t n_start = 0;
            
            for(size_t k = 0; k < niters; k++) {
              if(k > 1 && n_start < k - 1) n_start++;
              size_t nextNStart = n_start;
              
              for(size_t n = n_start; n < tend/dt; n++) {
                
                CLOG(INFO, "Parareal") << "n: " << n << " nstart: " << n_start;
                
                if(k == 0 || n > n_start) {
                  do_coarse(n > 0 ? u[n-1] : nullptr, coarseState, n*dt);
                  if(k < niters - 1) do_fine(n > 0 ? u[n-1] : nullptr, fineState, n*dt);
                }
                
                if(k > 0) {
                  diff->copy(u[n]);
                  u[n]->copy(uFine[n]);
                  
                  if(n > n_start) {
                    delta->copy(coarseState);
                    delta->saxpy(-1.0, uCoarse[n]);
                    // apply delta correction
                    u[n]->saxpy(1.0, delta);
                    CLOG(INFO, "Parareal") << "delta Norm: " << delta->norm0();
                  }
                  
                  diff->saxpy(-1.0,u[n]);
                  res = diff->norm0();
                  if(res < abs_res_tol) nextNStart++;
                }
                else {
                  u[n]->copy(coarseState);
                }
                uCoarse[n]->copy(coarseState);
                uFine[n]->copy(fineState);
                echo_error(n, k, res);
              } // loop over time slices
              n_start = nextNStart;
            } // loop over parareal iterations
           
            CLOG(INFO, "Parareal")  << "Time Measurement: " << double(clock() - timeMeasure)/CLOCKS_PER_SEC;
          } // function run_parareal
      }; // class Parareal
    } // namespace parareal
  } // namespace examples
} // namespace pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv,
               pfasst::examples::parareal::Parareal<>::init_opts,
               pfasst::examples::parareal::Parareal<>::init_logs);
  
  const double  dt     = pfasst::config::get_value<double>("dt", 0.01);
  const double  tend = pfasst::config::get_value<double>("tend", 0.1);
  const size_t  npariters = pfasst::config::get_value<size_t>("num_par_iter", size_t(tend/dt)+1);
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
  
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.run_parareal(tend, npariters, nfineiters, ncrseiters, 
                        dt, nnodes, nnodesCoarse, quad_type,
                        ndofs_fine, ndofs_coarse, abs_res_tol, rel_res_tol);
  
  fftw_cleanup();
}
#endif
