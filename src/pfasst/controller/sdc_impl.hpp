#include "pfasst/controller/sdc.hpp"
#include <iostream>

namespace pfasst
{
  template<typename time>
  void SDC<time>::run()
  {
    auto sweeper = this->get_level(0);

    std::cout << "StartTime: " << this->get_time() << " EndTime " << this->get_end_time() << "\n";
    for (; this->get_time() < this->get_end_time(); this->advance_time()) {
      bool initial = this->get_step() == 0;
      for (this->set_iteration(0);
           this->get_iteration() < this->get_max_iterations();
           this->advance_iteration()) {
        std::cout << "time: " << this->get_time() << " iteration: " << this->get_iteration() << "\n";
        bool predict = this->get_iteration() == 0;
        if (predict) {
          sweeper->predict(initial);
          sweeper->post_predict();
        } else {
          sweeper->sweep();
          sweeper->post_sweep();
        }
        if (sweeper->converged()) {
          break;
        }
      }
      sweeper->post_step();
      sweeper->advance();
    }
  }
}  // ::pfasst
