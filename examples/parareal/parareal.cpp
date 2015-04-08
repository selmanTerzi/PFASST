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
      class Parareal {

        private:
          vector<shared_ptr<SDC<time>>>  g; // coarse-level controllers
          vector<shared_ptr<SDC<time>>> f; // fine-level controllers
          
          shared_ptr<SpectralTransfer1D<time>> transfer;
          
          shared_ptr<Encapsulation<double>> uStart;
          
          vector<shared_ptr<Encapsulation<double>>> u; // vector with numerical solution
          vector<shared_ptr<Encapsulation<double>>> uExact; // vector with exact solution
          shared_ptr<Encapsulation<double>> err;
          error_map errors;
          
          shared_ptr<SDC<time>> getSDC(const size_t nnodes, const size_t ndofs, 
                                       const double dt, const double startTime, const bool fine) {
            auto sdc = make_shared<SDC<time>>();
            
            auto quad = quadrature::quadrature_factory(nnodes, pfasst::quadrature::QuadratureType::GaussLobatto);

            auto factory = make_shared<VectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            sweeper->set_residual_tolerances(1e-8, 0.0);
            sweeper->set_quadrature(quad);
            sweeper->set_factory(factory);

            sdc->add_level(sweeper);

            // 1 SDC controller for each time step
            sdc->set_duration(startTime, startTime+dt, dt, fine ? 10 : 1);
            sdc->setup();
            
            auto q0 = sweeper->get_start_state();
            sweeper->exact(q0, 0.0);
            
            return sdc;
          }
        
          void init(const size_t nsteps, const double dt, 
                    const size_t nnodes, const size_t ndofs) {
            
            size_t crsendofs = ndofs/2;
            
            auto factory = make_shared<VectorFactory<double>>(ndofs);
            auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
            
            g.clear();
            f.clear();
            for(size_t i = 0 ; i < nsteps; i++) {
//               g.push_back(getSDC(nnodes/3, crsendofs, dt, i*dt, false));
              f.push_back(getSDC(nnodes, ndofs, dt, i*dt, true));
              g.push_back(getSDC(nnodes, ndofs, dt, i*dt, true));
              
              u.push_back(factory->create(solution));
              uExact.push_back(factory->create(solution));
              sweeper->exact(uExact[i],i*dt);
            }
            err = factory->create(solution);
            uStart = factory->create(solution);
            transfer = make_shared<SpectralTransfer1D<>>();
          }
          
        public:
        
        
          static void init_opts()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts();
            pfasst::config::options::add_option<size_t>("Parareal", "num_par_iter", "num Parareal iterations");
          }
          
          static void init_logs()
          {
            pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_logs();
            pfasst::log::add_custom_logger("Parareal");
          }
        
          void run_parareal(const size_t nsteps, const size_t niters,
                            const double dt, const size_t nnodes, const size_t ndofs)
          {
            
            init(nsteps, dt, nnodes, ndofs);
            
            shared_ptr<SDC<time>> coarse;
            shared_ptr<SDC<time>> fine;
            
            shared_ptr<Encapsulation<double>> coarseState(nullptr);
            shared_ptr<Encapsulation<double>> deltaCourseState(nullptr);
            
            for(size_t k = 0; k < niters; k++) {
              for(size_t n = 0; n < nsteps; n++) {
                if(k > n + 1) continue;
                
                coarse = g[n];
                fine = f[n];
                
                std::cout << "k: " << k << " n: " << n << "\n";
                
                EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
                EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
                  
                if(k > 0) {
                  // calculate fine sweep and delta correction
                  CVLOG(1, "Parareal") << fineSweeper.get_start_state() << ": " << vector<double>(as_vector<double, time>(fineSweeper.get_start_state()));
                  CVLOG(1, "Parareal") << coarseSweeper.get_start_state() << ": " << vector<double>(as_vector<double, time>(coarseSweeper.get_start_state()));
                  if(n > 0) transfer->interpolate(fineSweeper.get_start_state(), coarseSweeper.get_start_state());
                  CVLOG(1, "Parareal") << fineSweeper.get_start_state() << ": " << vector<double>(as_vector<double, time>(fineSweeper.get_start_state()));
                  
                  u[n]->copy(coarseSweeper.get_start_state());
                  err->copy(uExact[n]);
                  err->saxpy(-1.0,u[n]);
                  std::cout << "Error: " << err->norm0() << "\n";
                  errors.insert(vtype(ktype(n, k-1), err->norm0()));
                  
                  std::cout << "Fine-Sweep" << std::endl;
                  fine->set_step(0);
                  fine->run();
                  
                  CLOG(INFO, "Parareal") << "Calculating delta";
                  CVLOG(4, "Parareal") << "verbose";
                  
                  if(!deltaCourseState) deltaCourseState = coarseSweeper.get_factory()->create(solution);
                  transfer->restrict(deltaCourseState, fineSweeper.get_end_state());
                  deltaCourseState->saxpy(-1.0, coarseSweeper.get_end_state());
                  
                  std::cout << "deltaCourseState Norm: " << deltaCourseState->norm0() << "\n";
                }
                
                if(coarseState) {
                  coarseSweeper.get_start_state()->copy(coarseState);
                  std::cout << "Coarse-StartState Copy, Norm: " << coarseState->norm0() << "\n";
                }
                else {
                  coarseState = coarseSweeper.get_factory()->create(solution);
                }
                
                std::cout << "Coarse-StartState Norm: " << coarseSweeper.get_start_state()->norm0() << "\n";
                
                std::cout << "Coarse-Sweep" << "\n"; 
                coarse->set_step(0);
                coarse->run();
                
                coarseState->copy(coarseSweeper.get_end_state());
                
                if(k > 0) {
                  // apply delta correction
                  coarseState->saxpy(1.0, deltaCourseState);
                  if(n == k-1) uStart->copy(coarseState);
                }
              } // loop over time steps
              if(k > 0) {
                coarseState->copy(uStart);
              }
              else {
                coarseState.reset();
              }
            } // loop over parareal iterations
          } // function run_parareal
      }; // Class Parareal
    } // namespace parareal
  } // namespace examples
} // namespace pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv,
               pfasst::examples::parareal::Parareal<>::init_opts,
               pfasst::examples::parareal::Parareal<>::init_logs);
  const size_t  nsteps = 8;
  const double  dt     = 0.01;
  const size_t  niters = 5;
  const size_t  nnodes = 9;
  const size_t  ndofs  = 64;
  
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.run_parareal(nsteps,niters,dt,nnodes,ndofs);
}
#endif
