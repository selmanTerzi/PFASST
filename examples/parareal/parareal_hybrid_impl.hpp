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
        shared_ptr<Encapsulation<time>> startState;
        
        auto finedelta = factory_fine->create(solution);
        auto crsedelta = factory_crse->create(solution);
        double res; // residual (difference of last and current iteration of end_state)
        
        bool prec_done = false; // boolean for checking if the precedessor is done
        bool done = false; // boolean for breaking next iteration if converged
        bool recvStartValue = false; // boolean for determining if a new startValue must be received
        
        CLOG(INFO, "Parareal") << "tend: " << this->get_end_time() 
                               << " dt: " << this->get_time_step()
                               << " num_iter: " << this->get_max_iterations();
        
        size_t nblocks = this->get_end_time()/this->get_time_step()/this->commSize+1;
        
        shared_ptr<ISweeper<>> fineSweeper = this->get_finest();
        shared_ptr<ISweeper<>> coarseSweeper = this->get_coarsest();
        
        MPI_Barrier(MPI_COMM_WORLD);
        double timeMeasure = MPI_Wtime();
        
        for(size_t nblock = 0; nblock < nblocks; nblock++) { // loop over time blocks
          this->set_step(commSize * nblock + commRank);
          
          CLOG(INFO, "Parareal") << "Time: " << this->get_time();
          if(this->get_time() >= this->get_end_time()) break;
          
          bool hasSuccessor = commRank < commSize - 1 &&
                              this->get_time() + this->get_time_step() <= this->get_end_time();
                              
          CLOG(INFO, "Parareal") << "hasSuccessor: " << hasSuccessor;
          for(this->set_iteration(0);
              this->get_iteration() < this->get_max_iterations() && !done; 
              this->advance_iteration()) { // loop over parareal iterations
            
            size_t k = this->get_iteration();
          
            bool predict = k == 0;
            
            if(!predict) {
              if(k == 1) {
                transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, true);
              }
              else if(recvStartValue) {
                transferFunc->interpolate(fineEncap->get_state(0), coarseEncap->get_state(0));
                fineEncap->reevaluate(true);
              }
              
              finedelta->copy(fineEncap->get_end_state());
              
              CLOG(INFO, "Parareal") << "Fine Sweep";
              fineEncap->sweep();
              fineEncap->post_sweep();
              
              // calc residium for break condition
              finedelta->saxpy(-1.0, fineEncap->get_end_state());
              res = finedelta->norm0();
              done = res < abs_res_tol;
              CLOG(INFO, "Parareal") << "Residual: " << res;
              if(done) CLOG(INFO, "Parareal") << "Done!";
              
              coarseEncap->save(false);
              transferFunc->PolyInterpMixin<time>::restrict(coarseSweeper, fineSweeper, true);
            }
            
            
            recvStartValue = !prec_done && (commRank > 0 || (nblock > 0 && predict));
            // get new initial value for coarse sweep
            if(recvStartValue) {
              int t = tag(k, nblock, commRank);
              coarseEncap->recv(comm, t, true);
              
              if(!predict && commRank > 0)
              {
                comm->status->recv(t);
                prec_done = comm->status->get_converged(commRank - 1);
              }
            }
            
            if(predict || hasSuccessor) {
              do_coarse(coarseState, predict);
            }
            
            // send new initial value to next processor
            if(hasSuccessor) {
              int t = tag(k, nblock,commRank+1);
              sendCorrection(crsedelta, coarseState, t);
              
              if(!predict) {
                comm->status->set_converged(done);
                comm->status->send(t);
              }
            }
          } // loop over parareal iterations
            
          comm->status->clear();
          prec_done = false;
          done = false;
          
          if(nblock < nblocks - 1 && !hasSuccessor) {
            transferFunc->restrict(coarseState, fineEncap->get_end_state());
            coarseState->send(comm, tag(0, nblock+1, 0), true);
          }
        } // loop over time blocks
          
        timeMeasure = MPI_Wtime() - timeMeasure;
        CLOG(INFO, "Parareal") << "time Measurement: " << timeMeasure;
      }
      
      template<typename time>
      void HybridParareal<time>::do_coarse(shared_ptr<Encapsulation<time>> end_state, 
                                           bool predict) 
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
        
        end_state->copy(coarseEncap->get_end_state());
      }
      
      template<typename time>
      void HybridParareal<time>::sendCorrection(shared_ptr<Encapsulation<time>> crsedelta,
                                                shared_ptr<Encapsulation<time>> coarseState, int tag)
      {
        // Calculate parareal-correction
        crsedelta->copy(coarseState);
        crsedelta->saxpy(-1.0, coarseEncap->get_saved_state(coarseEncap->get_nodes().size()-1));
        transferFunc->restrict(coarseState, fineEncap->get_end_state());
        coarseState->saxpy(1.0, crsedelta);
        coarseState->send(comm, tag, true);
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
