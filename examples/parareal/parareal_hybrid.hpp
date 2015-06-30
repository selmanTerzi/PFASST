#ifndef _PFASST__CONTROLLER__PARAREAL_HYBRID_HPP_
#define _PFASST__CONTROLLER__PARAREAL_HYBRID_HPP_

#include "pfasst/controller/interface.hpp"
#include "pfasst/encap/encap_sweeper.hpp"
#include "pfasst/encap/poly_interp.hpp"
#include "spectral_transfer_1d.hpp"

using namespace pfasst::encap;

namespace pfasst
{
  namespace examples 
  {
    namespace parareal
    {
      /**
      * Hybrid Parareal/SDC controller.
      *
      * @tparam time time precision
      *
      * @see Controller on how to set up the controller and feed it with sweepers to do the actual
      *   integration.
      *
      * @ingroup Controllers
      */
      template<typename time = time_precision>
      class HybridParareal
        : public Controller<time>
      {
        // MPI-variables
        ICommunicator* comm;
        int commSize;
        int commRank;
              
        shared_ptr<EncapFactory<>> factory_fine; // fine level factory
        shared_ptr<EncapFactory<>> factory_crse; // crse level factory
        
        EncapSweeper<time>* coarseEncap; // fine level sweeper
        EncapSweeper<time>* fineEncap; // coarse level sweeper 
        
        shared_ptr<SpectralTransfer1D<>> transferFunc; // transfer function
        
        double abs_res_tol; // residual tolerance for break condition
        
        private:
          virtual int tag(const size_t k, const size_t j, int commRank);
          
          void do_coarse(shared_ptr<Encapsulation<>> start_state,
                         shared_ptr<Encapsulation<>> end_state, 
                         bool predict);
          
          void recvState(shared_ptr<Encapsulation<double>> state, const size_t k, 
                         const size_t j, bool* prec_done);
          
        public:
          /**
          * Run hybrid parareal.
          */
          virtual void run();
          
          virtual void setup(double abs_res_tol, size_t ndofsfine, size_t ndofscoarse, 
                             shared_ptr<SpectralTransfer1D<>> transferFunc);
          
          virtual void set_comm(ICommunicator* comm);
      };
    }  // ::parareal
  }  // ::examples
}  // ::pfasst

#include "../examples/parareal/parareal_hybrid_impl.hpp"

#endif  // _PFASST__CONTROLLER__PARAREAL_HYBRID_HPP_
