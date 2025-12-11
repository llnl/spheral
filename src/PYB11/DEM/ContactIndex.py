#-------------------------------------------------------------------------------
# ContactIndex structure
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11module("SpheralCompiledModules.SpheralDEM")
class ContactIndex:
    def pyinit(self):
      "Constructor"
      return

    # Comparisons
    #def __eq__(self):
    #   return
    # def __lt__(self):
    #   return

    # String representation
    @PYB11implementation("""
[](const ContactIndex& self) {
  std::string result = "(" + 
                       std::to_string(self.storeNodeList) + " " +
                       std::to_string(self.storeNode) + " " +
                       std::to_string(self.storeContact) + " " +
                       std::to_string(self.pairNodeList) + " " +
                       std::to_string(self.pairNode) + " " +
                       std::to_string(self.solidBoundary) + ")";
  return result;
}""")
    def __repr__(self):
        return

    # Attributes
    storeNodeList = PYB11readwrite()
    storeNode = PYB11readwrite()
    storeContact = PYB11readwrite()
    pairNodeList = PYB11readwrite()
    pairNode = PYB11readwrite()
    solidBoundary = PYB11readwrite()

