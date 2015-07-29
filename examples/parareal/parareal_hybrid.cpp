#include <memory>
#include <vector>
using namespace std;

#include <mpi.h>
#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/mpi_communicator.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"
#include "parareal_hybrid.hpp"

namespace pfasst 
{
  namespace examples 
  {
    namespace parareal 
    {
      template<typename time = time_precision>
      shared_ptr<AdvectionDiffusionSweeper<time>> getSweeper(const size_t nnodes, 
                                       const pfasst::quadrature::QuadratureType quad_type, 
                                       const size_t ndofs, const double abs_res_tol, 
                                       const double rel_res_tol)
      {
        auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
        auto factory = make_shared<VectorFactory<double>>(ndofs);
        auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
        
        sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);
        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        
        return sweeper;
      }
      
      static void init_opts()
      {
        pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_opts();
        pfasst::config::options::add_option<size_t>("Parareal", "spatial_dofs_coarse", "Number of spatial degrees of freedom at coarse level");
        pfasst::config::options::add_option<size_t>("Parareal", "num_nodes_coarse", "Number of collocation nodes for coarse sweeper");
      }
      
      static void init_logs()
      {
        pfasst::examples::parareal::AdvectionDiffusionSweeper<>::init_logs();
        pfasst::log::add_custom_logger("Parareal");
      }
      
      void run_hybrid_parareal()
      {
        MPICommunicator comm(MPI_COMM_WORLD);
        HybridParareal<> para;
        
        const size_t num_nodes_fine = config::get_value<size_t>("num_nodes",5);
        const size_t num_nodes_coarse = config::get_value<size_t>("num_nodes_coarse",3);
        const size_t spatial_dofs = config::get_value<size_t>("spatial_dofs",128);
        const size_t spatial_dofs_coarse = config::get_value<size_t>("spatial_dofs_coarse",64);
        const double abs_res_tol = config::get_value<double>("abs_res_tol",1e-14);
        const double rel_res_tol = config::get_value<double>("rel_res_tol",1e-14);
        auto quadType = pfasst::quadrature::QuadratureType::GaussLobatto;
        
        auto transfer = make_shared<SpectralTransfer1D<>>();
        
        auto fineSweeper = getSweeper(num_nodes_fine, quadType, spatial_dofs, abs_res_tol, rel_res_tol);
        auto coarseSweeper = getSweeper(num_nodes_coarse, quadType, spatial_dofs_coarse, abs_res_tol, rel_res_tol);
        
        para.set_comm(&comm);
        para.add_level(fineSweeper);
        para.add_level(coarseSweeper);
        para.setup(abs_res_tol, spatial_dofs, spatial_dofs_coarse, transfer);
        para.set_options();
        
        fineSweeper->exact(fineSweeper->get_start_state(), 0.0);
        
        double timeMeasure = MPI_Wtime();
        para.run();
        timeMeasure = MPI_Wtime() - timeMeasure;
        CLOG(INFO, "Advec") << "time Measurement: " << timeMeasure;
      }
    }  // ::parareal
  }  // ::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  
  pfasst::init(argc, argv,
               pfasst::examples::parareal::init_opts,
               pfasst::examples::parareal::init_logs);
  
  pfasst::examples::parareal::run_hybrid_parareal();
  
  fftw_cleanup();
  MPI_Finalize();
}
#endif
