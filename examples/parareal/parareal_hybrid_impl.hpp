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
        
        CLOG(INFO, "Parareal") << "tend: " << this->get_end_time() 
                               << " dt: " << this->get_time_step()
                               << " num_iter: " << this->get_max_iterations();
        
        size_t numTiters = this->get_end_time()/this->get_time_step()/this->commSize;
        
        shared_ptr<ISweeper<>> fineSweeper = this->get_finest();
        shared_ptr<ISweeper<>> coarseSweeper = this->get_coarsest();
        
        MPI_Barrier(MPI_COMM_WORLD);
        double timeMeasure = MPI_Wtime();
        
        for(size_t j = 0; j < numTiters; j++) { // loop over time blocks
          
          size_t nglobal = commSize * j + commRank;
          if(this->get_time() > this->get_end_time()) break;
          
          this->set_step(nglobal);
          bool hasSuccessor = commRank < commSize - 1 &&
                              this->get_time() + this->get_time_step() <= this->get_end_time();
          
          for(this->set_iteration(0); 
              this->get_iteration() < this->get_max_iterations() && !done; 
              this->advance_iteration()) { // loop over parareal iterations
            
            size_t k = this->get_iteration();
          
            bool predict = k == 0;
            bool recvStartValue = commRank > 0 || (j > 0 && predict);
            
            if(!predict) {
              if(k == 1) {
                transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, true);
              }
              else if(commRank > 0) {
                transferFunc->interpolate(fineEncap->get_state(0), startState);
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
            
            // get new initial value for coarse sweep
            if(!prec_done && recvStartValue) {
              if(!startState) {
                startState = factory_crse->create(solution);
              }
              recvState(startState, k, j, &prec_done);
            }
            
            if(predict || hasSuccessor) {
              do_coarse(startState, coarseState, predict);
            }
            
            // send new initial value to next processor
            if(hasSuccessor) {
              int t = tag(k,j,commRank+1);
              
              // Calculate parareal-correction
              crsedelta->copy(coarseState);
              crsedelta->saxpy(-1.0, coarseEncap->get_saved_state(coarseEncap->get_nodes().size()-1));
              transferFunc->restrict(coarseState, fineEncap->get_end_state());
              coarseState->saxpy(1.0, crsedelta);
              coarseState->send(comm, t, true);
              
              if(!predict) {
                comm->status->set_converged(done);
                comm->status->send(t);
              }
            }
          } // loop over parareal iterations
            
          comm->status->clear();
          prec_done = false;
          done = false;
          
          if(j < numTiters - 1 && !hasSuccessor) {
            transferFunc->restrict(coarseState, fineEncap->get_end_state());
            coarseState->send(comm, tag(0,j+1,0), true);
          }
        } // loop over time blocks
          
        timeMeasure = MPI_Wtime() - timeMeasure;
        CLOG(INFO, "Parareal") << "time Measurement: " << timeMeasure;
      }
      
      template<typename time>
      void HybridParareal<time>::do_coarse(shared_ptr<Encapsulation<time>> start_state,
                                           shared_ptr<Encapsulation<time>> end_state, 
                                           bool predict) 
      {
        if(start_state) {
          if(predict) {
            coarseEncap->get_start_state()->copy(start_state);
          } 
          else {
            coarseEncap->get_state(0)->copy(start_state);
            coarseEncap->reevaluate(true);
          }
        }
        
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
      void HybridParareal<time>::recvState(shared_ptr<Encapsulation<time>> state, 
                                           const size_t k, const size_t j, bool* prec_done)
      {
        int t = tag(k, j, commRank);
        
        state->recv(comm, t, true);
        
        if(k > 0 && commRank > 0)
        {
          comm->status->recv(t);
          *prec_done = comm->status->get_converged(commRank - 1);
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
