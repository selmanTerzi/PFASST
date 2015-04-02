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
        
        shared_ptr<SDC<time>> getSDC(const size_t nnodes, const size_t ndofs, const double dt,
                                     const double startTime, const quadrature::QuadratureType quad_type) {
          auto sdc = make_shared<SDC<time>>();
          
          auto quad = quadrature::quadrature_factory(nnodes, quad_type);

          auto factory = make_shared<VectorFactory<double>>(ndofs);
          auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

          sweeper->set_quadrature(quad);
          sweeper->set_factory(factory);

          sdc->add_level(sweeper);

          // 1 SDC controller for each time step
          sdc->set_duration(startTime, startTime+dt, dt, 2);
          sdc->setup();
          
          auto q0 = sweeper->get_start_state();
          sweeper->exact(q0, 0.0);
          
          return sdc;
        }
        
        void init(const size_t nsteps, const double dt, const size_t nnodes, 
                  const size_t ndofs, const quadrature::QuadratureType quad_type) {
          g.clear();
          f.clear();
          for(size_t i = 0 ; i < nsteps; i++) {
            g.push_back(getSDC((nnodes-1)/3, ndofs/2, dt, i*dt, quad_type));
            f.push_back(getSDC(nnodes, ndofs, dt, i*dt, quad_type));
          }
          transfer = make_shared<SpectralTransfer1D<>>();
        }
        
      public:
        
        void run_parareal(const size_t nsteps, const size_t niters,
                          const double dt, const size_t nnodes, const size_t ndofs, 
                          const quadrature::QuadratureType quad_type) {
          
          init(nsteps, dt, nnodes, ndofs, quad_type);
          
          shared_ptr<SDC<time>> coarse;
          shared_ptr<SDC<time>> fine;
          
          shared_ptr<Encapsulation<double>> coarseState(nullptr);
          shared_ptr<Encapsulation<double>> fineState(nullptr);
          shared_ptr<Encapsulation<double>> deltaFineState(nullptr);
          
          for(size_t i = 0; i < niters; i++) {
            for(size_t j = 0; j < nsteps; j++) {
              if(i > j + 1) continue;
              
              coarse = *(g.begin()+j);
              fine = *(f.begin()+j);
              
              std::cout << "i: " << i << " j: " << j << "\n";
              
              EncapSweeper<time>& coarseSweeper = as_encap_sweeper<time>(coarse->get_finest());
              EncapSweeper<time>& fineSweeper = as_encap_sweeper<time>(fine->get_finest());
              
              if(i > 0) {
                if(j > 0) transfer->interpolate(fineSweeper.get_start_state(), coarseSweeper.get_start_state());
                fine->set_iteration(0);
                fine->set_step(0);
                fine->run();
                deltaFineState = fineSweeper.get_factory()->create(solution);
                deltaFineState->copy(fineSweeper.get_end_state());
                fineState = fineSweeper.get_factory()->create(solution);
                transfer->interpolate(fineState, coarseSweeper.get_end_state());
                deltaFineState->saxpy(-1.0, fineState);
              }
              
              if(coarseState) {
                coarseSweeper.get_start_state()->copy(coarseState);
              }
              else {
                coarseState = coarseSweeper.get_factory()->create(solution);
              }
              
              
              coarse->set_iteration(0);
              coarse->set_step(0);
              coarse->run();
              coarseState->copy(coarseSweeper.get_end_state());
              
              if(deltaFineState) {
                std::cout << "deltaFineState\n";
                transfer->interpolate(fineSweeper.get_end_state(), coarseState);
                fineSweeper.get_end_state()->saxpy(1.0, deltaFineState);
                transfer->restrict(coarseSweeper.get_end_state(),fineSweeper.get_end_state());
              }
            } // loop over time steps
            if(i > 0) {
              coarseState = as_encap_sweeper<time>((*(g.begin()+i-1))->get_finest()).get_end_state();
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
  const size_t  nsteps = 5;
  const pfasst::time_precision  dt     = 0.01;
  const size_t  niters = 3;
  const size_t  nnodes = 8;
  const size_t  ndofs  = 128;
  const pfasst::quadrature::QuadratureType quad_type = pfasst::quadrature::QuadratureType::GaussLegendre;
  
  pfasst::examples::parareal::Parareal<> parareal;
  parareal.run_parareal(nsteps,niters,dt,nnodes,ndofs,quad_type);
}
#endif
