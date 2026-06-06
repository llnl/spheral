#-------------------------------------------------------------------------------
# Adiak utilities
#-------------------------------------------------------------------------------
from PYB11Generator import *

# adiak::init() is called automatically when this module is loaded
# adiak::fini() is called automatically when this module is destroyed

@PYB11cppname("adiak::fini")
def adiak_fini():
    "Finish Adiak"
    return "void"

@PYB11cppname("adiak::collect_all")
def adiak_collect_all():
    "Add some default Adiak metadata"
    return "void"

adiak_categories = PYB11enum(("unset", "all", "general", "performance", "control"),
                             doc="Enum of Adiak categories")

@PYB11cppname("adiak::value")
@PYB11template("ValueType")
def adiak_value(name = "std::string",
                value = "%(ValueType)s",
                category = ("int",
                            "adiak_categories::general"),
                subcategory = ("std::string", '""')):
    "Set a single value in Adiak with a given name"
    return "bool"

@PYB11cppname("adiak::value")
@PYB11pyname("adiak_value")
@PYB11template("ValueType")
def adiak_value2(name = "std::string",
                 value = "%(ValueType)s",
                 value2 = "%(ValueType)s",
                 category = ("int",
                             "adiak_categories::general"),
                 subcategory = ("std::string", '""')):
    "Set a pair of values in Adiak with a given name"
    return "bool"
