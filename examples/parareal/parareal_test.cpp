/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <cstdlib>
#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>
#include <pfasst/encap/vector.hpp>
#include "parareal.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      error_map run_parareal(double abs_residual_tol, double rel_residual_tol=0.0)
      {
        PARAREAL<> parareal;

        auto const nnodes = config::get_value<size_t>("num_nodes", 3);
        auto const ndofs  = config::get_value<size_t>("spatial_dofs", 64);
        auto const quad_type = \
          config::get_value<quadrature::QuadratureType>("nodes_type", quadrature::QuadratureType::GaussLegendre);

        auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
        auto factory = make_shared<encap::VectorFactory<double>>(ndofs);
        auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        sweeper->set_residual_tolerances(abs_residual_tol, rel_residual_tol);

        parareal.add_level(sweeper);
        parareal.set_duration(0.0, 4*0.01, 0.01, 4);
        parareal.set_options();
        parareal.setup();

        auto q0 = sweeper->get_start_state();
        sweeper->exact(q0, 0.0);

        parareal.run();

        fftw_cleanup();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  pfasst::examples::advection_diffusion::run_parareal(0.0);
}
#endif
