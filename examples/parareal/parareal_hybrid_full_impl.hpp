#include "parareal_hybrid_full.hpp"

namespace pfasst
{
  namespace examples 
  {
    namespace parareal
    {
      template<typename time>
      void FullHybridParareal<time>::run()
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
        
        shared_ptr<ISweeper<time>> fineSweeper = this->get_finest();
        shared_ptr<ISweeper<time>> coarseSweeper = this->get_coarsest();
        
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
          
            auto u = factory_fine->create(solution);
            auto err = factory_fine->create(solution);
            auto advecSweeper = dynamic_pointer_cast<AdvectionDiffusionSweeper<>>(fineSweeper);
            advecSweeper->exact(u, this->get_time()+this->get_time_step());
          
            if(!predict) {
              err->copy(u);
              err->saxpy(-1.0, fineEncap->get_state(fineEncap->get_nodes().size()-1));
              CLOG(INFO, "Parareal") << "Fine Error before interpolate: " << err->norm0();
              
              CLOG(INFO, "Parareal") << "Interpolate";
              // compute parareal correction per interpolation
              if(k == 1) {
                transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, false);
              }
              else {
                transferFunc->interpolateDiff(fineSweeper, coarseSweeper);
              }
              
              err->copy(u);
              err->saxpy(-1.0, fineEncap->get_state(fineEncap->get_nodes().size()-1));
              CLOG(INFO, "Parareal") << "Fine Error after interpolate: " << err->norm0();
              
              CLOG(INFO, "Parareal") << "Fine Sweep";
              fineEncap->sweep();
              fineEncap->post_sweep();
              done = firstRank ? fineEncap->converged() : fineEncap->converged() && prec_done;
              if(done) CLOG(INFO, "Parareal") << "Done!";
              
              CLOG(INFO, "Parareal") << "Restrict";
              transferFunc->PolyInterpMixin<time>::restrict(coarseSweeper, fineSweeper, true);
            }
            
            bool recvStartValue = !prec_done && (!firstRank || (nblock > 0 && predict));
            if(recvStartValue) {
              CLOG(INFO, "Parareal") << "Receive";
              
              int t = tag(k, nblock, commRank);
              
              if(predict)
              {
                if(firstRank) {
                  CLOG(INFO, "Parareal") << "recv Fine Initial state";
                  fineEncap->recv(comm, t, true);
                }
                else {
                  CLOG(INFO, "Parareal") << "recv Coarse Initial state";
                  coarseEncap->recv(comm, t, true);
                }
              }
              else {
                CLOG(INFO, "Parareal") << "recv Coarse update state";
                coarseEncap->recv(comm, t, true);
                CLOG(INFO, "Parareal") << "Reevaluate coarse initial";
                coarseEncap->reevaluate(true);
                comm->status->recv(t);
                prec_done = comm->status->get_converged(commRank - 1);
              }
            }
            
            if(predict) {
              if(firstRank) {
                transferFunc->restrict_initial(coarseSweeper, fineSweeper);
                coarseEncap->get_state(0)->copy(coarseEncap->get_start_state());
              }
              else {
                transferFunc->interpolate_initial(fineSweeper, coarseSweeper);
              }
              fineEncap->get_state(0)->copy(fineEncap->get_start_state());
              fineEncap->spread();
              coarseEncap->spread();
              coarseEncap->save(false);
            }
            
            do_coarse(predict);
            
            // send new initial value to next processor
            if(!lastRank && hasSuccessor) {
              int t = tag(k, nblock,commRank+1);
              CLOG(INFO, "Parareal") << "Send coarse end_state";
              coarseEncap->send(comm, t, true);
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
            CLOG(INFO, "Parareal") << "Send Fine end_state to next block";
            fineEncap->send(comm, tag(0, nblock+1, 0), true);
          }
        } // loop over time blocks
      }
      
      template<typename time>
      void FullHybridParareal<time>::do_coarse(bool predict)
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
      int FullHybridParareal<time>::tag(const size_t k, const size_t j, int commRank) 
      {
        return k * 10000 + commSize * j + commRank;
      }
     
      template<typename time>
      void FullHybridParareal<time>::set_comm(ICommunicator* comm)
      {
        this->comm = comm;
        this->commRank = comm->rank();
        this->commSize = comm->size();
      }
      
      template<typename time>
      void FullHybridParareal<time>::setup(double abs_res_tol, size_t ndofsfine, size_t ndofscoarse, 
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
