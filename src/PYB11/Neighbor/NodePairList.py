from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodePairList
#-------------------------------------------------------------------------------
class NodePairList:
  def pyinit(self):
    "Default Constructor"

  def pyinit_vals(self,
                  vals = "const std::vector<NodePairIdxType>&"):
    "Constructs a NodePairList using a vector of NodePairIdxTypes"
    return

  def clear(self):
    "Clears all data from Node pair List."
    return "void"

  def size(self):
    "Returns the number of Node Pairs in the list"
    return "size_t"

  @PYB11implementation("[](const NodePairList& self) { return self.size(); }")
  def __len__(self):
    return

  @PYB11cppname("operator[]")
  @PYB11returnpolicy("reference_internal")
  @PYB11implementation('[](NodePairList& self, int i) { const int n = self.size(); return &self[(i %% n + n) %% n]; }')
  def __getitem__(self):
      return

  @PYB11implementation("[](const NodePairList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
  def __iter__(self):
    "Python iteration through a NodePairList."

  @PYB11implementation("""[](const NodePairList& self, const NodePairIdxType val) { for (const auto& ele: self) {
                                                                                      if (ele == val) return true;
                                                                                    }
                                                                                    return false;
                                                                                  }""")
  def __contains__(self,
                   val = "NodePairIdxType"):
    return "bool"

  @PYB11implementation("""[](const NodePairList& self, const py::tuple val0) { if (val0.size() != 4) throw py::value_error("require a tuple of 4 integers");
                                                                               NodePairIdxType val(val0[0].cast<int>(), val0[1].cast<int>(), val0[2].cast<int>(), val0[3].cast<int>());
                                                                               for (const auto& ele: self) {
                                                                                 if (ele == val) return true;
                                                                               }
                                                                               return false;
                                                                             }""")
  def __contains__(self,
                   val = "py::tuple"):
    return "bool"

