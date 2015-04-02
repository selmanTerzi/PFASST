#include <iostream>
#include <vector>

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
        
        shared_ptr<SDC<time>> getSDC(const size_t nnodes, const size_t ndofs, 
                                     const double dt, const double startTime) {
          auto sdc = make_shared<SDC<time>>();
          
          auto quad = quadrature::quadrature_factory(nnodes, pfasst::quadrature::QuadratureType::GaussLobatto);

          auto factory = make_shared<VectorFactory<double>>(ndofs);
          auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
          sweeper->set_residual_tolerances(1e-8, 0.0);

          sweeper->set_quadrature(quad);
          sweeper->set_factory(factory);

          sdc->add_level(sweeper);

          // 1 SDC controller for each time step
          sdc->set_duration(startTime, startTime+dt, dt, 10);
          sdc->setup();
          
          auto q0 = sweeper->get_start_state();
          sweeper->exact(q0, 0.0);
          
          return sdc;
        }
        
        void init(const size_t nsteps, const double dt, 
                  const size_t nnodes, const size_t ndofs) {
          g.clear();
          f.clear();
          for(size_t i = 0 ; i < nsteps; i++) {
//             g.push_back(getSDC((nnodes-1)/3, ndofs/2, dt, i*dt, quad_type));
            f.push_back(getSDC(nnodes, ndofs, dt, i*dt));
            g.push_back(getSDC(nnodes, ndofs, dt, i*dt));
          }
          transfer = make_shared<SpectralTransfer1D<>>();
        }
        
      public:
        
        void run_parareal(const size_t nsteps, const size_t niters,
                          const double dt, const size_t nnodes, const size_t ndofs) {
          
          init(nsteps, dt, nnodes, ndofs);
          
          shared_ptr<SDC<time>> coarse;
          shared_ptr<SDC<time>> fine;
          
          shared_ptr<Encapsulation<double>> coarseState(nullptr);
          shared_ptr<Encapsulation<double>> fineState(nullptr);
          shared_ptr<Encapsulation<double>> deltaFineState(nullptr);
          shared_ptr<Encapsulation<double>> debugState(nullptr);
          
          for(size_t k = 0; k < niters; k++) {
            for(size_t n = 0; n < nsteps; n++) {
              if(k > n + 1) continue;
              
              coarse = *(g.begin()+n);
              fine = *(f.begin()+n);
              
              std::cout << "k: " << k << " n: " << n << "\n";
              
              EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
              EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
              
              // DEBUG:
              if(!debugState) debugState = coarseSweeper.get_factory()->create(solution);
              debugState->copy(coarseSweeper.get_start_state());
              debugState->saxpy(-1.0, coarseSweeper.get_state(0));
              std::cout << "Debug-Error Coarse: " << debugState->norm0() << "\n";
              debugState->copy(fineSweeper.get_start_state());
              debugState->saxpy(-1.0, fineSweeper.get_state(0));
              std::cout << "Debug-Error Fine: " << debugState->norm0() << "\n";
                
              if(k > 0) {
                if(n > 0) transfer->interpolate(fineSweeper.get_start_state(), coarseSweeper.get_start_state());
                
                std::cout << "Fine-Sweep" << "\n";
                
                fine->set_step(0);
                fine->run();
                
                if(!deltaFineState) deltaFineState = fineSweeper.get_factory()->create(solution);
                deltaFineState->copy(fineSweeper.get_end_state());
                
                if(!fineState) fineState = fineSweeper.get_factory()->create(solution);
                transfer->interpolate(fineState, coarseSweeper.get_end_state());
                
                deltaFineState->saxpy(-1.0, fineState);
              }
              
              if(coarseState) {
                coarseSweeper.get_start_state()->copy(coarseState);
              }
              else {
                coarseState = coarseSweeper.get_factory()->create(solution);
              }
              
              std::cout << "Coarse-Sweep" << "\n";
              coarse->set_step(0);
              coarse->run();
              coarseState->copy(coarseSweeper.get_end_state());
              
              if(deltaFineState) {
                transfer->interpolate(fineState, coarseState);
                fineState->saxpy(1.0, deltaFineState);
                transfer->restrict(coarseState,fineState);
                coarseSweeper.get_end_state()->copy(coarseState);
              }
            } // loop over time steps
            if(k > 0) {
              coarseState->copy(as_encap_sweeper<time>((*(g.begin()+k-1))->get_finest()).get_end_state());
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
               pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_logs);
  const size_t  nsteps = 4;
  const double  dt     = 0.01;
  const size_t  niters = 2;
  const size_t  nnodes = 3;
  const size_t  ndofs  = 64;
  
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.run_parareal(nsteps,niters,dt,nnodes,ndofs);
}
#endif
