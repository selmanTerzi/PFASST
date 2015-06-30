/*
 * Spectral (FFT) transfer routines.
 */

#ifndef _SPECTRAL_TRANSFER_1D_HPP_
#define _SPECTRAL_TRANSFER_1D_HPP_

#include <cassert>
#include <cstdlib>
#include <memory>
using namespace std;

#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/poly_interp.hpp>
#include <pfasst/encap/encap_sweeper.hpp>

#include "fft.hpp"



namespace pfasst
{
  namespace examples
  {
    namespace parareal
    {
      template<typename time = pfasst::time_precision>
      class SpectralTransfer1D
        : public encap::PolyInterpMixin<time>
      {
          typedef encap::Encapsulation<double> Encapsulation;

          FFT fft;

        public:
          void interpolate(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src) override
          { 
            auto& fine = encap::as_vector<double, time>(dst);
            auto& crse = encap::as_vector<double, time>(src);
            
            if(fine.size() == crse.size()) {
              dst->copy(src);
              return;
            }
            
            auto* crse_z = this->fft.forward(crse);
            auto* fine_z = this->fft.get_workspace(fine.size())->z;

            for (size_t i = 0; i < fine.size(); i++) {
              fine_z[i] = 0.0;
            }

            double c = 1.0 / crse.size();

            for (size_t i = 0; i < crse.size() / 2; i++) {
              fine_z[i] = c * crse_z[i];
            }

            for (size_t i = 1; i < crse.size() / 2; i++) {
              fine_z[fine.size() - crse.size() / 2 + i] = c * crse_z[crse.size() / 2 + i];
            }

            this->fft.backward(fine);
          }
          
          // interpolates the difference of coarse solution and last fine solution at coarse nodes
          // to the fine nodes
          void interpolateDiff(shared_ptr<ISweeper<time>> dst,
                               shared_ptr<const ISweeper<time>> src)
          {
            auto& fine = encap::as_encap_sweeper(dst);
            auto& crse = encap::as_encap_sweeper(src);

            if (this->tmat.rows() == 0) {
              this->tmat = pfasst::quadrature::compute_interp<time>(fine.get_nodes(), crse.get_nodes());
            }
            
            auto const crse_nodes = crse.get_nodes();
            auto const fine_nodes = fine.get_nodes();
            
            size_t nfine = fine_nodes.size();
            size_t ncrse = crse_nodes.size();
            
            auto crse_factory = crse.get_factory();
            auto fine_factory = fine.get_factory();

            vector<shared_ptr<Encapsulation>> diff(ncrse), fine_state(nfine);
            shared_ptr<Encapsulation> restricted = crse_factory->create(encap::EncapType::solution);
            shared_ptr<Encapsulation> diffCoarse = crse_factory->create(encap::EncapType::solution);
            
            for (size_t m = 0; m < nfine; m++) { fine_state[m] = fine.get_state(m); }
            
            int trat = (int(nfine) - 1) / (int(ncrse) - 1);
            
            for (size_t m = 0; m < ncrse; m++) {
              diffCoarse->copy(crse.get_state(m));
              
              if (crse_nodes[m] != fine_nodes[m * trat]) {
                throw NotImplementedYet("coarse nodes must be nested");
              }
              this->restrict(restricted, fine_state[m * trat]);
              diffCoarse->saxpy(-1.0, restricted);
              diff[m] = fine_factory->create(encap::EncapType::solution);
              this->interpolate(diff[m], diffCoarse);
            }
            
            // interpolate the difference of coarse solution and restricted fine solution at coarse nodes to fine nodes
            fine.get_state(0)->mat_apply(fine_state, 1.0, this->tmat, diff, false);
            
            fine.reevaluate();
          }

          void restrict(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(src);
            auto& crse = encap::as_vector<double, time>(dst);

            size_t xrat = fine.size() / crse.size();

            for (size_t i = 0; i < crse.size(); i++) {
              crse[i] = fine[xrat*i];
            }
          }
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#endif
