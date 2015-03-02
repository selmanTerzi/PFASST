#ifndef _PFASST_PARAREAL_HPP
#define _PFASST_PARAREAL_HPP

#include "pfasst/controller/interface.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace parareal
    {
      /**
       * Parareal controller
       */
      template<typename time = time_precision>
      class PARAREAL
        : public Controller<time>
      {
        protected:
          typedef typename pfasst::Controller<time>::LevelIter LevelIter;
          
        public :
          virtual void run();
      };
      
    }
  }
}

#include "parareal_impl.hpp"

#endif