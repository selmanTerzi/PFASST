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
        shared_ptr<Encapsulation<double>> startState;
        
        auto diff = factory_fine->create(solution);
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
          if(nglobal*this->get_time_step() > this->get_end_time()) break;
          
          for(this->set_iteration(0); 
              this->get_iteration() < this->get_max_iterations() && !done; 
              this->advance_iteration()) { // loop over parareal iterations
            
            size_t k = this->get_iteration();
          
            bool predict = k == 0;
            if(!predict) {
              if(k == 1) {
                CLOG(INFO, "Parareal") << "Interpolating coarse to fine";
                transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, true);
              }
              else {
                CLOG(INFO, "Parareal") << "Interpolating difference from coarse to fine";
                transferFunc->interpolateDiff(fineSweeper, coarseSweeper);
              }
              coarseEncap->save(false);
              diff->copy(fineEncap->get_end_state());
              CLOG(INFO, "Parareal") << "Fine Sweep";
              fineEncap->sweep();
              fineEncap->post_sweep();
              CLOG(INFO, "Parareal") << "Restricting from fine to coarse";
              transferFunc->PolyInterpMixin<time>::restrict(coarseSweeper, fineSweeper, true);
              coarseEncap->reevaluate();
            }
            
            // get new initial value for coarseEncap sweep
            if(!prec_done && (commRank > 0 || (predict && j > 0))) {
              if(!startState) startState = factory_crse->create(solution);
              recvState(startState, k, j, &prec_done);
            }
            
            do_coarse(startState, coarseState, predict);
            
            // calc residium for break condition
            if(!predict) {
              // calc the residium
              diff->saxpy(-1.0, fineEncap->get_end_state());
              res = diff->norm0();
              done = res < abs_res_tol;
              
              CLOG(INFO, "Parareal") << "Residual: " << res;
              CLOG(INFO, "Parareal") << "Done: " << done;
            }
            
            // send new initial value to next processor
            if(commRank < commSize - 1 && (nglobal+1)*this->get_time_step() <= this->get_end_time()) {
              int t = tag(k, j, commRank+1);
              CLOG(INFO, "Parareal") << "Send tag: " << t;
              if(!predict) transferFunc->restrict(coarseState, fineEncap->get_end_state());
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
          
          if(j < numTiters - 1 && commRank == commSize - 1) {
            CLOG(INFO, "Parareal") << "Send tag: " << tag(0, j+1, 0);
            transferFunc->restrict(coarseState, fineEncap->get_end_state());
            coarseState->send(comm, tag(0, j+1, 0), true); // Send to proc with rank 0
          }
        } // loop over time blocks
          
        timeMeasure = MPI_Wtime() - timeMeasure;
        CLOG(INFO, "Parareal") << "time Measurement: " << timeMeasure;
      }
      
      template<typename time>
      void HybridParareal<time>::do_coarse(shared_ptr<Encapsulation<>> start_state,
                                           shared_ptr<Encapsulation<>> end_state, 
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
            
        if(predict) CLOG(INFO, "Parareal") << "Predict";
        else CLOG(INFO, "Parareal") << "Coarse-Sweep";
        
        if(predict) {
          coarseEncap->predict(true);
          coarseEncap->post_predict();
        }
        else {
          coarseEncap->sweep();
          coarseEncap->post_sweep();
        } 
        
        end_state->copy(coarseEncap->get_end_state());
      }
      
      template<typename time>
      void HybridParareal<time>::recvState(shared_ptr<Encapsulation<double>> state, 
                                           const size_t k, const size_t j, bool* prec_done)
      {
        int t = tag(k, j, commRank);
        CLOG(INFO, "Parareal") << "Recv tag: " << t;
        
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
        
        this->coarseEncap = &encap::as_encap_sweeper<time>(this->get_coarsest());
        this->fineEncap = &encap::as_encap_sweeper<time>(this->get_finest());
        this->transferFunc = transferFunc;
        
        this->fineEncap->set_controller(this);
        this->coarseEncap->set_controller(this);
        this->fineEncap->setup(false);
        this->coarseEncap->setup(true);
        
        this->factory_fine = make_shared<VectorFactory<double>>(ndofsfine);
        this->factory_crse = make_shared<VectorFactory<double>>(ndofscoarse);
      }
    }  // ::parareal
  }  // ::examples
}  // ::pfasst
