from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodePairIdxType
#-------------------------------------------------------------------------------
class NodePairIdxType:
    def pyinit(self,
               i_n = "int",
               i_l = "int",
               j_n = "int",
               j_l = "int",
               f = ("double", "1.0")):
      "Constructor"
      return

    # Comparisons
    def __eq__(self):
      return
    def __lt__(self):
      return

    # String representation
    @PYB11implementation("""
[](const NodePairIdxType& self) {
  std::string result = "(" + 
                       std::to_string(self.i_list) + " " +
                       std::to_string(self.i_node) + " " +
                       std::to_string(self.j_list) + " " +
                       std::to_string(self.j_node) + ")";
  return result;
}""")
    def __repr__(self):
        return

    # Attributes
    i_node = PYB11readwrite()
    j_node = PYB11readwrite()
    i_list = PYB11readwrite()
    j_list = PYB11readwrite()
    f_couple = PYB11readwrite()
