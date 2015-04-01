
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
        list<shared_ptr<SDC<time>>>  g; // coarse-level controllers
        list<shared_ptr<SDC<time>>> f; // fine-level controllers
        
        const ITransfer<time> transfer;
        
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
          sdc->set_duration(startTime, startTime+dt, dt, 1);
          sdc->setup();
          
          auto q0 = sweeper->get_start_state();
          sweeper->exact(q0, 0.0);
          
          return sdc;
        }
        
        void init(const size_t nsteps, const double dt, const size_t nnodes, 
                  const size_t ndofs, const quadrature::QuadratureType quad_type) {
          g.clear();
          f.clear();
          for(int i= 0 ; i< nsteps; i++) {
            g.insert(getSDC((nnodes-1)/3, ndofs/2, dt, i*dt, quad_type));
            f.insert(getSDC(nnodes, ndofs, dt, i*dt, quad_type));
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
          
          for(int i = 0; i < niters; i++) {
            for(int j = 0; j < nsteps; j++) {
              if(i > j + 1) continue;
              
              coarse = *(g.begin()+j);
              fine = *(f.begin()+j);
              
              EncapSweeper<time> coarseSweeper = as_encap_sweeper<time>(coarse->finest());
              EncapSweeper<time> fineSweeper = as_encap_sweeper<time>(fine->finest());
              
              if(i > 0) {
                transfer.interpolate(fineSweeper.get_start_state(), coarseSweeper.get_start_state());
                fine.set_step(0);
                fine.run();
                deltaFineState->copy(fineSweeper.get_end_state());
                transfer.interpolate(fineState, coarseSweeper.get_end_state());
                deltaFineState->saxpy(-1.0, fineState);
              }
              
              if(coarseState) {
                coarseSweeper.get_start_state()->copy(coarseState);
              }
              
              coarse.set_step(0);
              coarse.run();
              coarseState->copy(coarseSweeper.get_end_state());
              
              if(deltaFineState) {
                transfer.interpolate(fineSweeper.get_end_state(), coarseState);
                fineSweeper.get_end_state()->saxpy(1.0, deltaFineState);
                transfer.restrict(coarseSweeper.get_end_state(),fineSweeper.get_end_state());
              }
            }
            if(i > 0) {
              coarseState = as_encap_sweeper<time>((*(g.begin()+i-1))->finest()).get_end_state();
            }
            else {
              coarseState.reset();
            }
          }
        }
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
  auto nsteps = 16;
  auto dt     = 0.01;
  auto niters = 4;
  auto const nnodes = 8;
  auto const ndofs  = 128;
  auto const quad_type = pfasst::quadrature::QuadratureType::GaussLegendre;
        
  auto parareal = make_shared<pfasst::examples::parareal::Parareal<pfasst::time_precision>>();
  parareal->run_parareal(nsteps,dt,nnodes,niters,ndofs,quad_type);
  }
#endif
