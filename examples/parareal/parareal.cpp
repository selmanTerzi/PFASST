#include <iostream>
#include <memory>
#include <vector>
using namespace std;

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/encap_sweeper.hpp>

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
          
          double dt;
          size_t nfineiters;
          size_t ncrseiters;
          
          shared_ptr<SDC<time>> getSDC(const size_t nnodes, const size_t ndofs)
          {
            auto sdc = make_shared<SDC<time>>();
            
            auto quad = quadrature::quadrature_factory(nnodes, pfasst::quadrature::QuadratureType::GaussLobatto);
            
            auto factory = make_shared<VectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            sweeper->set_residual_tolerances(1e-8, 0.0);
            sweeper->set_quadrature(quad);
            sweeper->set_factory(factory);
            
            sdc->add_level(sweeper);
            sdc->setup();
            
            auto q0 = sweeper->get_start_state();
            sweeper->exact(q0, 0.0);
            
            return sdc;
          }
        
          void init(const size_t nsteps, const double dt, const size_t nnodes, 
                    const size_t ndofs_fine, const size_t ndofs_coarse,
                    const size_t nfineiters, const size_t ncrseiters) 
          {
            
            factory = make_shared<VectorFactory<double>>(ndofs_fine);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs_fine);
            this->dt = dt;
            this->nfineiters = nfineiters;
            this->ncrseiters = ncrseiters;
            
            fine = getSDC(nnodes, ndofs_fine);
            coarse = getSDC(nnodes, ndofs_coarse);
//            coarse = getSDC(nnodes/3, ndofs_coarse);
            
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
            
//             CVLOG(1, "Parareal") << "End_state coarse" << vector<double>(as_vector<double,time>(coarseSweeper.get_end_state()));
            
            if(as_vector<double>(coarseSweeper.get_end_state()).size() == as_vector<double>(end_state).size()) {
              end_state->copy(coarseSweeper.get_end_state());
            }
            else {
              transfer->interpolate(end_state, coarseSweeper.get_end_state());
            }
            
//             CVLOG(1, "Parareal") << "End_state interpolated" << vector<double>(as_vector<double,time>(end_state));
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
          
        public:
        
        
          static void init_opts()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts();
            pfasst::config::options::add_option<size_t>("Parareal", "num_par_iter", "Number of Parareal iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "spatial_dofs_coarse", "Number of spatial degrees of freedom at coarse level");
            pfasst::config::options::add_option<size_t>("Parareal", "num_fine_iter", "Number of Fine Sweep iterations");
            pfasst::config::options::add_option<size_t>("Parareal", "num_crse_iter", "Number of Coarse Sweep iterations");
          }
          
          static void init_logs()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_logs();
            pfasst::log::add_custom_logger("Parareal");
          }
        
          void run_parareal(const size_t nsteps, const size_t niters,
                            const size_t nfineiters, const size_t ncrseiters,
                            const double dt, const size_t nnodes, 
                            const size_t ndofs_fine, const size_t ndofs_coarse)
          {
            
            init(nsteps, dt, nnodes, ndofs_fine, ndofs_coarse, nfineiters, ncrseiters);
            
            shared_ptr<Encapsulation<double>> coarseState = factory->create(solution);
            shared_ptr<Encapsulation<double>> fineState = factory->create(solution);
            shared_ptr<Encapsulation<double>> delta = factory->create(solution);
            
            for(size_t k = 0; k < niters; k++) {
              for(size_t n = 0; n < nsteps; n++) {
                if(k > n + 1) continue;
                
                CLOG(INFO, "Parareal") << "k: " << k << " n: " << n;
                
                if(n >= k) do_coarse(n > 0 ? u[n-1] : nullptr, coarseState, n*dt);
                if(n >= k) do_fine(n > 0 ? u[n-1] : nullptr, fineState, n*dt);
                
                if(k > 0) {
                  u[n]->copy(uFine[n]);
                  
                  if(n >= k) {
                    CLOG(INFO, "Parareal") << "Calculating delta";
                    
                    delta->copy(coarseState);
                    delta->saxpy(-1.0, uCoarse[n]);
                    
                    CLOG(INFO, "Parareal") << "delta Norm: " << delta->norm0();
                    
                    // apply delta correction
                    u[n]->saxpy(1.0, delta);
                  }
                }
                else {
                  u[n]->copy(coarseState);
                }
                uCoarse[n]->copy(coarseState);
                uFine[n]->copy(fineState);
                
                err->copy(uExact[n]);
                err->saxpy(-1.0, u[n]);
                CLOG(INFO, "Parareal") << "Error: " << err->norm0();
                errors.insert(vtype(ktype(n, k-1), err->norm0()));
              } // loop over time slices
            } // loop over parareal iterations
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
  const size_t  nsteps = 5;
  const double  dt     = 0.01;
  const size_t  npariters = pfasst::config::get_value<size_t>("num_par_iter", 5);
  const size_t  nfineiters = pfasst::config::get_value<size_t>("num_fine_iter", 10);
  const size_t  ncrseiters = pfasst::config::get_value<size_t>("num_crse_iter", 5);
  const size_t  nnodes = pfasst::config::get_value<size_t>("num_nodes", 9);
  const size_t  ndofs_fine  = pfasst::config::get_value<size_t>("spatial_dofs", 64);
  const size_t  ndofs_coarse  = pfasst::config::get_value<size_t>("spatial_dofs_coarse", 64);
  
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.run_parareal(nsteps, npariters, nfineiters, ncrseiters, dt, nnodes, ndofs_fine, ndofs_coarse);
}
#endif
