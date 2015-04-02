#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  inline static void init(int argc, char** argv,
                          std::function<void()> opts = nullptr,
                          std::function<void()> logs = nullptr)
  {
    if (opts) {
      opts();
    }
    config::init();
    log::start_log(argc, argv);
    if (logs) {
      logs();
    }
    config::read_commandline(argc, argv);
  }
} // ::pfasst


/**
 * @defgroup Controllers
 *   Controllers represent the main entry point of PFASST++ as they ensemble the central algorithmic
 *   logic of PFASST and related algorithms.
 *
 * @defgroup Internals
 *   Entities listed in this module are ment to be for internal use only and should not be used
 *   outside of PFASST++ itself.
 *
 * @defgroup Utilities
 *   General utility functions not directly related to PFASST++ but also handy in user code.
 *
 * @defgroup Examples
 *   A few different examples demonstrating the use and functionality of PFASST++.
 */

#endif
