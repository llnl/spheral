#-------------------------------------------------------------------------------
# TimerMgr class
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11singleton
class TimerMgr:

    "Singleton wrapper for CaliperManager. Accesses a C++ singleton object."

    @PYB11static
    def add(self, config_str = "std::string"):
        "Add a Caliper configuration"
        return "void"

    @PYB11static
    def load(self, config_json = "std::string"):
        "Load a json file containing Caliper configurations"
        return "void"

    @PYB11static
    def default_start(self, testname = "std::string"):
        "Set the spot Caliper configuration and start the manager"
        return "void"

    @PYB11static
    def start(self):
        "Start the Caliper configuration manager"
        return "void"

    @PYB11static
    def stop(self):
        "Stop the Caliper configuration manager"
        return "void"

    @PYB11static
    def fini(self):
        "Flush the Caliper configuration manager"
        return "void"

    @PYB11static
    def get_config(self):
        "Return the current Caliper configuration"
        return "std::string"

    @PYB11static
    def get_filename(self):
        "Return current Caliper filename, if set"
        return "std::string"

    @PYB11static
    def timers_usable(self):
        "Return whether the code has been compiled with timers turned on"
        return "bool"
    @PYB11static
    def timer_begin(self, region_name = "std::string"):
        "Start custom region Caliper timer, must have corresponding timer_end call"
        return "void"

    @PYB11static
    def timer_end(self, region_name = "std::string"):
        "End custom region Caliper timer"
        return "void"

    @PYB11static
    def is_started(self):
        "Check if ConfigManager has been started"
        return "bool"
