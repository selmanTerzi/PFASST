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
        // variables for the residium calculation via the difference of the current and last iterations fine end state
        auto finedelta = factory_fine->create(solution);
        double res; // residual (difference of last and current iteration of end_state)
        
        bool firstRank = commRank == 0;
        bool lastRank = commRank == commSize -1;
        bool prec_done = false; // boolean for checking if the precedessor is done
        bool done = false; // boolean for breaking next iteration if converged
        bool recvStartValue = false; // boolean for determining if a new startValue must be received
        double timeMeasure = 0.0; // variable for timings
        double tInterpolRestrict = 0.0; // variable for timing of interpolation and restriction operations
        double tCommunication = 0.0; // variable for timing of communication
        
        double div = this->get_end_time()/this->get_time_step();
        if(div - size_t(div) > 0) {
          CLOG(ERROR, "Controller") << "invalid time step: dt must be a divisor of tend";
          throw ValueError("invalid time step: dt must be a divisor of tend");
        }
        size_t nblocks = div/commSize - size_t(div)/commSize > 0 ? size_t(div)/commSize + 1 : size_t(div)/commSize;
        
        shared_ptr<ISweeper<time>> fineSweeper = this->get_finest();
        shared_ptr<ISweeper<time>> coarseSweeper = this->get_coarsest();
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(size_t nblock = 0; nblock < nblocks; nblock++) { // loop over time blocks
          this->set_step(commSize * nblock + commRank);
          
          CVLOG(2, "Parareal") << "Time: " << this->get_time();
          if(this->get_time() >= this->get_end_time()) break;
          
          bool initial = firstRank && nblock == 0;
          bool hasSuccessor = this->get_time() + this->get_time_step() < this->get_end_time();
          
          for(this->set_iteration(0);
              this->get_iteration() < this->get_max_iterations() && !done; 
              this->advance_iteration()) { // loop over parareal iterations
            
            size_t k = this->get_iteration();
            bool predict = k == 0;
          
            if(!predict) {
              CVLOG(2, "Parareal") << "Interpolate";
              // compute parareal correction per interpolation
              timeMeasure = MPI_Wtime();
              transferFunc->PolyInterpMixin<time>::interpolate(fineSweeper, coarseSweeper, true);
              tInterpolRestrict += MPI_Wtime() - timeMeasure;
              
              // set finedelta to last fine end state for the calculation of the difference
              if(this->diffResidual) finedelta->copy(fineEncap->get_end_state());
              
              CLOG(INFO, "Parareal") << "Fine Sweep";
              fineEncap->sweep();
              fineEncap->post_sweep();
              
              if(this->diffResidual) {
                // calc residium for break condition (difference of current and last fine end state):
                finedelta->saxpy(-1.0, fineEncap->get_end_state());
                res = finedelta->norm0();
                done = res < abs_res_tol;
                CLOG(INFO, "Parareal") << "DiffResidual: " << res;
              }
              else {
                CVLOG(2, "Parareal") << "fine converged:" << fineEncap->converged();
                done = fineEncap->converged();
              }
              if(done) CVLOG(2, "Parareal") << "Done!";
              
              CVLOG(2, "Parareal") << "Restrict";
              timeMeasure = MPI_Wtime();
              transferFunc->PolyInterpMixin<time>::restrict(coarseSweeper, fineSweeper, true);
              tInterpolRestrict += MPI_Wtime() - timeMeasure;
              
              coarseEncap->save(false);
            }
            
            recvStartValue = !prec_done && (!firstRank || (nblock > 0 && predict));
            // get new initial value for coarse sweep
            if(recvStartValue) {
              int t = tag(k, nblock, commRank);
              
              CVLOG(2, "Parareal") << "recv Coarse Initial state";
              timeMeasure = MPI_Wtime();
              coarseEncap->recv(comm, t, true);
              
              if(!predict) {
                CVLOG(2, "Parareal") << "Reevaluate coarse initial";
                coarseEncap->reevaluate(true);
                comm->status->recv(t);
                prec_done = comm->status->get_converged(commRank - 1);
              }
              tCommunication += MPI_Wtime() - timeMeasure;
            }
            
            if(predict) {
              if(initial) {
                CVLOG(2, "Parareal") << "restrict initial state";
                timeMeasure = MPI_Wtime();
                transferFunc->restrict_initial(coarseSweeper, fineSweeper);
                tInterpolRestrict += MPI_Wtime() - timeMeasure;
              }
              else {
                CVLOG(2, "Parareal") << "interpolate initial state";
                timeMeasure = MPI_Wtime();
                transferFunc->interpolate_initial(fineSweeper, coarseSweeper);
                tInterpolRestrict += MPI_Wtime() - timeMeasure;
                
                if(fineEncap->get_quadrature()->left_is_node()) {
                  fineEncap->get_state(0)->copy(fineEncap->get_start_state());
                }
              }  
              CVLOG(2, "Parareal") << "fine spread";
              fineEncap->spread();
              CVLOG(2, "Parareal") << "coarse spread";
              coarseEncap->spread();
              coarseEncap->save(false);
            }
            
            do_coarse(predict);
            
            // send new initial value to next processor
            if(!lastRank && hasSuccessor) {
              int t = tag(k, nblock,commRank+1);
              CVLOG(2, "Parareal") << "Send coarse end_state";
              timeMeasure = MPI_Wtime();
              coarseEncap->send(comm, t, true);
              if(!predict) {
                comm->status->set_converged(done);
                comm->status->send(t);
              }
              tCommunication += MPI_Wtime() - timeMeasure;
            }
          } // loop over parareal iterations
            
          comm->status->clear();
          prec_done = false;
          done = false;
          
          if(lastRank && nblock < nblocks - 1 && hasSuccessor) {
            CVLOG(2, "Parareal") << "Send restricted fine end_state to next block";
            timeMeasure = MPI_Wtime();
            transferFunc->restrict(coarseEncap->get_end_state(), fineEncap->get_end_state());
            tInterpolRestrict += MPI_Wtime() - timeMeasure;
            
            
            timeMeasure = MPI_Wtime();
            coarseEncap->send(comm, tag(0, nblock+1, 0), true);
            tCommunication += MPI_Wtime() - timeMeasure;
          }
        } // loop over time blocks
        CLOG(INFO, "Advec") << "time Measurement Interpolation: " << tInterpolRestrict;
        CLOG(INFO, "Advec") << "time Measurement Communication: " << tCommunication;
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
                                           shared_ptr<SpectralTransfer1D<>> transferFunc, bool diffResidual)
      {
        this->abs_res_tol = abs_res_tol;
        this->transferFunc = transferFunc;
        
        this->factory_fine = make_shared<VectorFactory<time>>(ndofsfine);
        this->factory_crse = make_shared<VectorFactory<time>>(ndofscoarse);
        
        this->coarseEncap = &encap::as_encap_sweeper<time>(this->get_coarsest());
        this->fineEncap = &encap::as_encap_sweeper<time>(this->get_finest());
        
        this->diffResidual = diffResidual;
        
        fineEncap->set_controller(this);
        coarseEncap->set_controller(this);
        fineEncap->setup(false);
        coarseEncap->setup(true);
      }
    }  // ::parareal
  }  // ::examples
}  // ::pfasst
