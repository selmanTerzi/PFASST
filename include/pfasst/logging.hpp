#ifndef _PFASST__LOGGING_HPP_
#define _PFASST__LOGGING_HPP_

#include "site_config.hpp"

#include <string>
using namespace std;

#include "config.hpp"

// enable easy logging of STL containers
#define _ELPP_STL_LOGGING
// disable creation of default log file
#define _ELPP_NO_DEFAULT_LOG_FILE
// enable passing `--logging-flags` via command line
#define _ELPP_LOGGING_FLAGS_FROM_ARG

#ifndef NDEBUG
  #define _ELPP_DEBUG_ASSERT_FAILURE
  #define _ELPP_STACKTRACE_ON_CRASH
#endif

#include <easylogging++.h>

#ifndef PFASST_LOGGER_INITIALIZED
  // initialize easyloggingpp
  // FIXME: this might already be called by code using PFASST++
  _INITIALIZE_EASYLOGGINGPP
  /**
   * guard symbol to ensure easylogging++ is only initialized once
   *
   * When this symbol is defined, it is expected that `_INITIALIZE_EASYLOGGINGPP` has been called once to initialize
   * easylogging++.
   *
   * \note In case the executable using PFASST++ is also using easylogging++ as its logging library and initializing it
   *     prior to including the PFASST++ headers, please define this symbol prior including any of the PFASST++ headers.
   */
  #define PFASST_LOGGER_INITIALIZED
#endif

#ifndef PFASST_LOGGER_DEFAULT_GLOBAL_FORMAT
  //! format for the default global logger
  #define PFASST_LOGGER_DEFAULT_GLOBAL_FORMAT "[%level] %msg"
#endif

#ifndef PFASST_LOGGER_DEFAULT_GLOBAL_TOFILE
  //! whether to log to file by default
  #define PFASST_LOGGER_DEFAULT_GLOBAL_TOFILE "false"
#endif
#ifndef PFASST_LOGGER_DEFAULT_GLOBAL_TOSTDOUT
  //! whether to log to stdout by default
  #define PFASST_LOGGER_DEFAULT_GLOBAL_TOSTDOUT "true"
#endif

#ifndef PFASST_LOGGER_DEFAULT_DEBUG_FORMAT
  //! format for the default global debug logger
  #define PFASST_LOGGER_DEFAULT_DEBUG_FORMAT "[%datetime{%Y-%M-%d %H:%m:%s,%g}] %level - %fbase:%line - %msg"
#endif

namespace pfasst
{
  namespace log
  {

    /**
     * sets default configuration for default loggers
     */
    static void load_default_config()
    {
      el::Configurations defaultConf;
      defaultConf.setToDefault();
      defaultConf.set(el::Level::Global, el::ConfigurationType::Format, PFASST_LOGGER_DEFAULT_GLOBAL_FORMAT);
      defaultConf.set(el::Level::Global, el::ConfigurationType::ToFile, PFASST_LOGGER_DEFAULT_GLOBAL_TOFILE);
      defaultConf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput, PFASST_LOGGER_DEFAULT_GLOBAL_TOSTDOUT);
      defaultConf.set(el::Level::Debug, el::ConfigurationType::Format, PFASST_LOGGER_DEFAULT_DEBUG_FORMAT);
      el::Loggers::reconfigureLogger("default", defaultConf);
    }

    /**
     * sets some default flags for easylogging++
     *
     * Current defaults are:
     *
     * - NewLineForContainer
     * - LogDetailedCrashReason
     * - ColoredTerminalOutput
     * - MultiLoggerSupport
     * - HierarchicalLogging
     * - AutoSpacing
     *
     * \see https://github.com/easylogging/easyloggingpp#logging-flags
     */
    static void set_logging_flags()
    {
      el::Loggers::addFlag(el::LoggingFlag::NewLineForContainer);
      el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
      el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
      el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
      el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
      el::Loggers::addFlag(el::LoggingFlag::AutoSpacing);
    }

    /**
     * starts easylogging++ with given arguments and loads configuration
     *
     * Usually, you want to pass the command line arguments given to `main()` in here and let easylogging++ figure
     * out what it needs.
     *
     * \param[in] argc number of command line arguments
     * \param[in] argv command line arguments
     */
    static void start_log(int argc, char** argv)
    {
      _START_EASYLOGGINGPP(argc, argv);
      load_default_config();
      set_logging_flags();
    }
  }  // ::pfasst::log
}  // ::pfasst

#endif  // _PFASST__LOGGING_HPP_