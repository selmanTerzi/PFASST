#include "parareal.hpp"

using pfasst::examples::parareal::PARAREAL;

namespace pfasst
{
  namespace examples
  {
    namespace parareal
    {
      /**
       * Parareal controller
       */
      template<typename time>
      void PARAREAL<time>::run()
      {
        auto coarseSweeper = this->get_coarsest();
        auto fineSweeper = this->get_finest();
        auto trns = this->get_transfer((size_t)0);
       
        for (; this->get_time() < this->get_end_time(); this->advance_time()) {
          bool initial = this->get_step() == 0;
          for (this->set_iteration(0);
               this->get_iteration() < this->get_max_iteration();
               this->advance_iteration()) {
            bool predict = this->get_iteration() == 0;
            if (predict) {
              coarseSweeper->predict(initial);
              coarseSweeper->post_predict();
              trns->interpolate(fineSweeper, coarseSweeper);
              fineSweeper->sweep();
              fineSweeper->post_sweep();
            } else {
              trns-> fas(this->get_time_step(), coarseSweeper, fineSweeper);
              coarseSweeper->sweep();
              coarseSweeper->post_sweep();
              trns->interpolate(fineSweeper, coarseSweeper);
              fineSweeper->sweep();
              fineSweeper->post_sweep();
              
            }
            coarseSweeper->save();
            fineSweeper->save();
            
            if ( fineSweeper->converged()) {
              break;
            }
            
            coarseSweeper->post_step();
            coarseSweeper->advance();
          }
        }
      }
    }
  }
}