from PYB11Generator import *

@PYB11template("DataType")
@PYB11holder("std::shared_ptr")
class IncrementalStatistic:
    def pyinit(self,
               meanGuess = ("const %(DataType)s", 0),
               name = ("const std::string", "\"IncrementalStatistic\""),
               printAdd = ("bool", "false")):
        "Allows tracking of incrementally added statistics"
        
    def add(self,
            data = "const %(DataType)s"):
        return "void"

    @PYB11const
    def mean(self):
        return "%(DataType)s"
    @PYB11const
    def variance(self):
        return "%(DataType)s"
    @PYB11const
    def total(self):
        return "%(DataType)s"
    @PYB11const
    def min(self):
        return "%(DataType)s"
    @PYB11const
    def max(self):
        return "%(DataType)s"
    @PYB11const
    def meanGuess(self):
        return "%(DataType)s"
    @PYB11const
    def numPoints(self):
        return "%(DataType)s"
    @PYB11const
    def shiftedSum(self):
        return "%(DataType)s"
    @PYB11const
    def shiftedSum2(self):
        return "%(DataType)s"
    @PYB11cppname("print")
    @PYB11const
    def print1(self):
        return "void"
