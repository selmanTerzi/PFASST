/**
 * Advection-Diffusion with vanilla SDC.
 *
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/vanilla_sdc.cpp
 * @since v0.1.0
 */
#include <cstdlib>
#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Advection/diffusion example using an encapsulated IMEX sweeper.
       *
       * This example uses a vanilla SDC sweeper.
       *
       * @ingroup AdvectionDiffusion
       */
      error_map run_vanilla_sdc(double abs_residual_tol, double rel_residual_tol=0.0)
      {
        SDC<> sdc;

        auto const dt = config::get_value<double>("dt", 0.01);
        auto const nnodes = config::get_value<size_t>("num_nodes", 3);
        auto const ndofs  = config::get_value<size_t>("spatial_dofs", 64);
        auto const tend = config::get_value<double>("tend", 0.04);
        auto const quad_type = \
          config::get_value<quadrature::QuadratureType>("nodes_type", quadrature::QuadratureType::GaussLobatto);
        
        const size_t niters = pfasst::config::get_value<double>("niters", 10);
        
        auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
        auto factory = make_shared<encap::VectorFactory<double>>(ndofs);
        auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        sweeper->set_residual_tolerances(abs_residual_tol, rel_residual_tol);

        sdc.add_level(sweeper);
        sdc.set_duration(0.0, tend, dt, niters);
        sdc.set_options();
        sdc.setup();

        auto q0 = sweeper->get_start_state();
        sweeper->exact(q0, 0.0);
        
        clock_t timeMeasure = clock();
        sdc.run();
        CLOG(INFO, "Advec")  << "Time Measurement: " << double(clock() - timeMeasure)/CLOCKS_PER_SEC;
        
        fftw_cleanup();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 1e-10);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);
  
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  pfasst::examples::advection_diffusion::run_vanilla_sdc(abs_res_tol, rel_res_tol);
}
#endif
