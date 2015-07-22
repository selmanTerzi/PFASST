#include "../examples/parareal/parareal_hybrid.hpp"


namespace pfasst
{
  namespace examples 
  {
    namespace parareal
    {
      template<typename time>
      void HybridParareal<time>::run()
      {
        if(!comm || commSize == 1) {
          CLOG(ERROR, "Parareal") << "Number of processors must be greater than 1";
          MPI_Finalize();
          exit(-1);
        }
        auto coarseState = factory_crse->create(solution);
        
        bool firstRank = commRank == 0;
        bool lastRank = commRank == commSize -1;
        bool prec_done = false; // boolean for checking if the precedessor is done
        bool done = false; // boolean for breaking next iteration if converged
        
        double div = this->get_end_time()/this->get_time_step();
        size_t nblocks = div - size_t(div) > 0 ? size_t(div) + 1 : size_t(div);
        nblocks = double(nblocks)/commSize - nblocks/commSize > 0 ? nblocks/commSize + 1 : nblocks/commSize;
        
        shared_ptr<ISweeper<>> fineSweeper = this->get_finest();
        shared_ptr<ISweeper<>> coarseSweeper = this->get_coarsest();
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(size_t nblock = 0; nblock < nblocks; nblock++) { // loop over time blocks
          this->set_step(commSize * nblock + commRank);
          
          CLOG(INFO, "Parareal") << "Time: " << this->get_time();
          if(this->get_time() >= this->get_end_time()) break;
          
          bool hasSuccessor = this->get_time() + this->get_time_step() < this->get_end_time();
          
          for(this->set_iteration(0);
              this->get_iteration() < this->get_max_iterations() && !done; 
              this->advance_iteration()) { // loop over parareal iterations
            
            size_t k = this->get_iteration();
            bool predict = k == 0;
            
            if(!predict) {
              CLOG(INFO, "Parareal") << "Fine Sweep";
              fineEncap->sweep();
              fineEncap->post_sweep();
              done = firstRank ? fineEncap->converged() : fineEncap->converged() && prec_done;
            }
            
            bool recvStartValue = !prec_done && (!firstRank || (nblock > 0 && predict));
            if(recvStartValue) {
              int t = tag(k, nblock, commRank);
              coarseEncap->recv(comm, t, true);
              if(lastRank) {
                transferFunc->interpolate_initial(fineSweeper, coarseSweeper);
                if(fineEncap->get_quadrature()->left_is_node()) {
                  fineEncap->get_state(0)->copy(fineEncap->get_start_state());
                  fineEncap->reevaluate(true);
                }
              }
              if(!predict && !firstRank)
              {
                comm->status->recv(t);
                prec_done = comm->status->get_converged(commRank - 1);
              }
            }
            
            bool doCoarse = predict || (!firstRank && !lastRank && recvStartValue);
            if(doCoarse) {
              do_coarse(predict);
              CLOG(INFO, "Parareal") << "Interpolate";
              transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, true);
              coarseEncap->save(false);
            }
            
            // send new initial value to next processor
            if(!lastRank && hasSuccessor) {
              int t = tag(k, nblock,commRank+1);
              if(predict) {
                coarseEncap->send(comm, t, true);
              }
              else {
                transferFunc->restrict(coarseState, fineEncap->get_end_state());
                coarseState->send(comm, t, true);
              }
              if(!predict) {
                comm->status->set_converged(done);
                comm->status->send(t);
              }
            }
          } // loop over parareal iterations
            
          comm->status->clear();
          prec_done = false;
          done = false;
          
          if(lastRank && nblock < nblocks - 1 && hasSuccessor) {
            CLOG(INFO, "Parareal") << "Send Restricted Fine end_state to next block";
            transferFunc->restrict(coarseState, fineEncap->get_end_state());
            coarseState->send(comm, tag(0, nblock+1, 0), true);
          }
        } // loop over time blocks
      }
      
      template<typename time>
      void HybridParareal<time>::do_coarse(bool predict)
      {
        if(predict) {
          CLOG(INFO, "Parareal") << "Coarse Predict";
          coarseEncap->predict(true);
          coarseEncap->post_predict();
        }
        else {
          CLOG(INFO, "Parareal") << "Coarse Sweep";
          coarseEncap->sweep();
          coarseEncap->post_sweep();
        }
      }
      
      template<typename time>
      int HybridParareal<time>::tag(const size_t k, const size_t j, int commRank) 
      {
        return k * 10000 + commSize * j + commRank;
      }
     
      template<typename time>
      void HybridParareal<time>::set_comm(ICommunicator* comm)
      {
        this->comm = comm;
        this->commRank = comm->rank();
        this->commSize = comm->size();
      }
      
      template<typename time>
      void HybridParareal<time>::setup(double abs_res_tol, size_t ndofsfine, size_t ndofscoarse, 
                                       shared_ptr<SpectralTransfer1D<>> transferFunc)
      {
        this->abs_res_tol = abs_res_tol;
        this->transferFunc = transferFunc;
        
        this->factory_fine = make_shared<VectorFactory<time>>(ndofsfine);
        this->factory_crse = make_shared<VectorFactory<time>>(ndofscoarse);
        
        this->coarseEncap = &encap::as_encap_sweeper<time>(this->get_coarsest());
        this->fineEncap = &encap::as_encap_sweeper<time>(this->get_finest());
        
        fineEncap->set_controller(this);
        coarseEncap->set_controller(this);
        fineEncap->setup(false);
        coarseEncap->setup(true);
      }
    }  // ::parareal
  }  // ::examples
}  // ::pfasst
