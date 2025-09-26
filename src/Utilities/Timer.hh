//---------------------------------Spheral++----------------------------------//
// Timer macros and Caliper ConfigManger singleton
//
//----------------------------------------------------------------------------//
#ifndef __Spheral_Timer__
#define __Spheral_Timer__
#include "config.hh"
// If SPHERAL_ENABLE_TIMERS is defined then we want timer functionality
#ifdef SPHERAL_ENABLE_TIMERS
#include "caliper/cali.h"
#include "caliper/cali-manager.h"

#define TIME_FUNCTION CALI_CXX_MARK_FUNCTION
#define TIME_SCOPE(regionName) CALI_CXX_MARK_SCOPE(regionName)
#define TIME_BEGIN(regionName) CALI_MARK_BEGIN(regionName)
#define TIME_END(regionName) CALI_MARK_END(regionName)
// For larger code regions, use phase
#define TIME_PHASE_BEGIN(regionName) CALI_MARK_PHASE_BEGIN(regionName)
#define TIME_PHASE_END(regionName) CALI_MARK_PHASE_END(regionName)

#else // SPHERAL_ENABLE_TIMERS
// Stub TIME macros, when SPHERAL_ENABLE_TIMERS is off

#define TIME_FUNCTION
#define TIME_SCOPE(regionName)
#define TIME_BEGIN(regionName)
#define TIME_END(regionName)
#define TIME_PHASE_BEGIN(regionName)
#define TIME_PHASE_END(regionName)

#endif // SPHERAL_ENABLE_TIMERS
// Note: This class is initialized in
// SimulationControl/SpheralTimingParser.py
namespace Spheral {
class TimerMgr {
private:
  TimerMgr() = default;
  ~TimerMgr() { }
  TimerMgr(const TimerMgr&) = delete;
  TimerMgr& operator=(const TimerMgr&) = delete;
  bool started = false;
  std::string caliperFilename = "";
  std::string caliperConfig = "";
public:
  static TimerMgr& instance() {
    static TimerMgr theInstance;
    return theInstance;
  }
  static void timer_begin(std::string regionName) {
    TIME_BEGIN(regionName.c_str());
  }
  static void timer_end(std::string regionName) {
    TIME_END(regionName.c_str());
  }
  static bool is_started() {
    return instance().started;
  }
  static std::string get_config() {
    return instance().caliperConfig;
  }
  static std::string get_filename() {
    return instance().caliperFilename;
  }
#ifdef SPHERAL_ENABLE_TIMERS
private:
  cali::ConfigManager cali_mgr;
public:
  static bool timers_usable() {
    return true;
  }
  static void add(std::string config_str) {
    bool test = instance().cali_mgr.add(config_str.c_str());
    VERIFY2(test, instance().cali_mgr.error_msg());
    instance().caliperConfig += config_str;
  }
  static void load(std::string config_json) {
    instance().cali_mgr.load(config_json.c_str());
    VERIFY2(!instance().cali_mgr.error(), instance().cali_mgr.error_msg());
  }
  static void default_start(std::string testname) {
    if (!testname.empty()) {
      std::string default_config = "spot,output=" + testname;
      instance().caliperFilename = testname;
      add(default_config);
      start();
    } else if (Spheral::Process::getRank() == 0) {
      std::cout << "WARNING: Caliper test name is empty, "
                << "no Caliper configuration started" << std::endl;
    }
  }
  static void start() {
    instance().cali_mgr.start();
    instance().started = true;
  }
  static void stop() {
    instance().cali_mgr.stop();
  }
  static void fini() {
    instance().cali_mgr.flush();
    instance().started = false;
  }
#else
  static bool timers_usable() {
    return false;
  }
  static void load(std::string) {
  }
  static void default_start(std::string) {
  }
  static void add(std::string) {
  }
  static void start() {
  }
  static void stop() {
  }
  static void fini() {
  }
#endif
};
}
#endif
