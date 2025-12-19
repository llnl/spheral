"""
Provides wrappers for the Silo library.
"""

from PYB11Generator import *

PYB11includes = ['"Geometry/Dimension.hh"',
                 '"Geometry/GeomPlane.hh"',
                 '"SiloWrappers.hh"']

PYB11namespaces = ["silo",
                   "Spheral"]

PYB11opaque = ["std::vector<char>",
               "std::vector<unsigned>",
               "std::vector<uint64_t>",
               "std::vector<int>",
               "std::vector<float>",
               "std::vector<double>",
               "std::vector<std::string>",

               "std::vector<std::vector<char>>",
               "std::vector<std::vector<unsigned>>",
               "std::vector<std::vector<uint64_t>>",
               "std::vector<std::vector<int>>",
               "std::vector<std::vector<float>>",
               "std::vector<std::vector<double>>",
               "std::vector<std::vector<std::string>>"]

#-------------------------------------------------------------------------------
class DBfile:
    "Opaque object for silo file struct"

#-------------------------------------------------------------------------------
@PYB11cppname("DBoptlist_wrapper")
class DBoptlist:
    "The silo optlist collection."

    # Constructor
    def pyinit(self,
               maxopts = ("int", "10")):
        "Construct with the given number of avilable option slots."

    #...........................................................................
    # addOption
    # @PYB11template("T")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            %(T)s& value) {
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(&value));
    #                            }""")
    # def addOption(self, option="int", value="%(T)s&"):
    #     return "int"

    # @PYB11template("T")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            int option_size,
    #                            std::vector<%(T)s>& value) {
    #                                auto value_size = value.size();
    #                                auto result = DBAddOption(self.mOptlistPtr, option_size, static_cast<void*>(&value_size));
    #                                if (result != 0) return result;
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(&(value.front())));
    #                            }""")
    # def addOptionVec(self, option="int", value="std::vector<%(T)s>&"):
    #     return "int"

    # addOptionInt =       PYB11TemplateMethod(addOption,    template_parameters="int",    pyname="addOption")
    # addOptionDouble =    PYB11TemplateMethod(addOption,    template_parameters="double", pyname="addOption")
    # addOptionVecInt =    PYB11TemplateMethod(addOptionVec, template_parameters="int",    pyname="addOption")
    # addOptionVecDouble = PYB11TemplateMethod(addOptionVec, template_parameters="double", pyname="addOption")

    # @PYB11pyname("addOption")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            std::string& value) {
    #                                return DBAddOption(self.mOptlistPtr, option, const_cast<char*>(value.c_str()));
    #                            }""")
    # def addOptionString(self, option="int", value="std::string&"):
    #     return "int"

    # @PYB11pyname("addOption")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            int option_size,
    #                            std::vector<std::string>& value) {
    #                                auto value_size = value.size();
    #                                auto result = DBAddOption(self.mOptlistPtr, option_size, static_cast<void*>(&value_size));
    #                                if (result != 0) return result;
    #                                char** charArray = new char*[value.size()];
    #                                for (auto k = 0; k < value.size(); ++k) {
    #                                  charArray[k] = new char[value[k].size() + 1];
    #                                  strcpy(charArray[k], value[k].c_str());
    #                                }
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(charArray));
    #                            }""")
    # def addOptionVecString(self, option="int", value="std::vector<std::string>&"):
    #     return "int"

    for (T, Label) in (("int", "Int"),
                       ("double", "Double"),
                       ("std::string", "String")):
        exec('''
@PYB11pycppname("addOption")
def addOption%(Label)s(self,
                       option = "int",
                       value = "const %(T)s&"):
    return "int"

@PYB11pycppname("getOption")
def getOption%(Label)s(self,
                       option = "int"):
            return "%(T)s"

@PYB11pycppname("addOption")
def addOptionVec%(Label)s(self,
                          option = "int",
                          option_size = "int",
                          value = "const std::vector<%(T)s>&"):
    return "int"

@PYB11pycppname("getOption")
def getOptionVec%(Label)s(self,
                          option = "int",
                          option_size = "int"):
    return "std::vector<%(T)s>"
''' % {"T" : T,
       "Label" : Label})

#-------------------------------------------------------------------------------
@PYB11cppname("DBmrgtree_wrapper")
class DBmrgtree:
    "Another silo thingus."

    # Constructor
    def pyinit(self,
               mesh_type = ("int", "DB_POINTMESH"),
               info_bits = ("int", "0"),
               max_children = ("int", "1024"),
               optlist = ("silo::DBoptlist_wrapper", "silo::DBoptlist_wrapper(10)")):
        "Constructor"

    # Methods for defining properties
    @PYB11const
    @PYB11cppname("name")
    @PYB11ignore
    def getname(self):
        return "std::string"

    @PYB11cppname("name")
    @PYB11ignore
    def setname(self, val="std::string"):
        return "void"

    @PYB11const
    @PYB11cppname("src_mesh_name")
    @PYB11ignore
    def getsrc_mesh_name(self):
        return "std::string"

    @PYB11cppname("src_mesh_name")
    @PYB11ignore
    def setsrc_mesh_name(self, val="std::string"):
        return "void"

    @PYB11const
    @PYB11cppname("src_mesh_type")
    @PYB11ignore
    def getsrc_mesh_type(self):
        return "int"

    @PYB11cppname("src_mesh_type")
    @PYB11ignore
    def setsrc_mesh_type(self, val="int"):
        return "void"

    @PYB11const
    @PYB11cppname("type_info_bits")
    @PYB11ignore
    def gettype_info_bits(self):
        return "int"

    @PYB11cppname("type_info_bits")
    @PYB11ignore
    def settype_info_bits(self, val="int"):
        return "void"

    @PYB11const
    @PYB11cppname("num_nodes")
    @PYB11ignore
    def getnum_nodes(self):
        return "int"

    @PYB11cppname("num_nodes")
    @PYB11ignore
    def setnum_nodes(self, val="int"):
        return "void"

    # Properties
    name = property(getname, setname)
    src_mesh_name = property(getsrc_mesh_name, setsrc_mesh_name)
    src_mesh_type = property(getsrc_mesh_type, setsrc_mesh_type)
    type_info_bits = property(gettype_info_bits, settype_info_bits)
    num_nodes = property(getnum_nodes, setnum_nodes)

#-------------------------------------------------------------------------------
class DBdefvars:
    ndefs = PYB11readwrite(doc="number of definitions")
    names = PYB11property(getterraw="[](DBdefvars& self) { return copy2py(self.names, self.ndefs); }",
                          setterraw="[](DBdefvars& self, py::list vals) { copy2c(self.names, vals); }",
                          doc="[ndefs] derived variable names")
    types = PYB11property(getterraw="[](DBdefvars& self) { return copy2py(self.types, self.ndefs); }",
                          setterraw="[](DBdefvars& self, py::list vals) { copy2c(self.types, vals); }",
                          doc="[ndefs] derived variable types")
    defns = PYB11property(getterraw="[](DBdefvars& self) { return copy2py(self.defns, self.ndefs); }",
                          setterraw="[](DBdefvars& self, py::list vals) { copy2c(self.defns, vals); }",
                          doc="[ndefs] derived variable definitions")
    guihides = PYB11property(getterraw="[](DBdefvars& self) { return copy2py(self.names, self.ndefs); }",
                          setterraw="[](DBdefvars& self, py::list vals) { copy2c(self.guihides, vals); }",
                          doc="[ndefs] flags to hide from post-processor's GUI")

#-------------------------------------------------------------------------------
class DBpointmesh:
    id = PYB11readwrite(doc="Identifier for this object")
    block_no = PYB11readwrite(doc="Block number for this mesh")
    group_no = PYB11readwrite(doc="Block group number for this mesh")
    name = PYB11readwrite(doc="Name associated with this mesh")
    cycle = PYB11readwrite(doc="Problem cycle number")
    units = PYB11property(getterraw="[](DBpointmesh& self) { return copy2py(self.units, 3); }",
                          setterraw="[](DBpointmesh& self, py::list vals) { copy2c(self.units, vals); }",
                          doc="Units for each axis")
    labels = PYB11property(getterraw="[](DBpointmesh& self) { return copy2py(self.labels, 3); }",
                           setterraw="[](DBpointmesh& self, py::list vals) { copy2c(self.labels, vals); }",
                           doc="Labels for each axis")
    title = PYB11readwrite(doc="Title for curve")
    coords = PYB11property(getterraw = "[](DBpointmesh& self) { return copy2py(self.coords, 3, self.nels, self.datatype); }",
                           setterraw="[](DBpointmesh& self, py::list vals) { copy2c(self.coords, vals, 3, self.nels, self.datatype); }",
                           doc = "Coordinate values")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    datatype = PYB11readwrite(doc="Datatype for coords (float, double)")
    ndims = PYB11readwrite(doc="Number of computational dimensions")
    nels = PYB11readwrite(doc="Number of elements in mesh")
    origin = PYB11readwrite(doc="0' or '1'")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    mrgtree_name = PYB11readwrite(doc="optional name of assoc. mrgtree object")
    gnznodtype = PYB11readwrite(doc="datatype for global node/zone ids")
    ghost_node_labels = PYB11readwrite()

#-------------------------------------------------------------------------------
class DBmultimesh:
    "A silo multimesh"

    # Attributes
    id = PYB11readwrite(doc="Identifier for this object")
    nblocks = PYB11readwrite(doc="Number of blocks in mesh")
    ngroups = PYB11readwrite(doc="Number of block groups in mesh")
    meshids = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.meshids, self.nblocks); }",
                            setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.meshids, vals); }",
                            doc="Array of mesh-ids which comprise mesh")
    meshnames = PYB11property(getterraw = "[](DBmultimesh& self) { return copy2py(self.meshnames, self.nblocks); }",
                              setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.meshnames, vals); }",
                              doc="Array of mesh-names for meshids")
    meshtypes = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.meshtypes, self.nblocks); }",
                              setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.meshtypes, vals); }",
                              doc="Array of mesh-type indicators [nblocks]")
    dirids = PYB11readwrite(doc="Array of directory IDs which contain blk")
    blockorigin = PYB11readwrite(doc="Origin (0 or 1) of block numbers")
    grouporigin = PYB11readwrite(doc="Origin (0 or 1) of group numbers")
    extentssize = PYB11readwrite(doc="size of each extent tuple")
    extents = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.extents, 2*self.nblocks); }",
                            setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.extents, vals); }",
                            doc="min/max extents of coords of each block")
    zonecounts = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.zonecounts, self.nblocks); }",
                               setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.zonecounts, vals); }",
                               doc="array of zone counts for each block")
    has_external_zones = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.has_external_zones, self.nblocks); }",
                                       setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.has_external_zones, vals); }",
                                       doc="external flags for each block")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    lgroupings = PYB11readwrite(doc="size of groupings array")
    groupings = PYB11readwrite(doc="Array of mesh-ids, group-ids, and counts")
    groupnames = PYB11property(getterraw = "[](DBmultimesh& self) { return copy2py(self.groupnames, self.ngroups); }",
                               setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.groupnames, vals); }",
                               doc="Array of group-names for groupings")
    mrgtree_name = PYB11readwrite(doc="optional name of assoc. mrgtree object")
    tv_connectivity = PYB11readwrite()
    disjoint_mode = PYB11readwrite()
    topo_dim = PYB11readwrite(doc="Topological dimension; max of all blocks.")
    file_ns = PYB11readwrite(doc="namescheme for files (in lieu of meshnames)")
    block_ns = PYB11readwrite(doc="namescheme for block objects (in lieu of meshnames)")
    block_type = PYB11readwrite(doc="constant block type for all blocks (in lieu of meshtypes)")
    empty_list = PYB11property(getterraw="[](DBmultimesh& self) { return copy2py(self.empty_list, self.empty_cnt); }",
                               setterraw="[](DBmultimesh& self, py::list vals) { copy2c(self.empty_list, vals); }",
                               doc="list of empty block #'s (option for namescheme)")
    empty_cnt = PYB11readwrite(doc="size of empty list")
    repr_block_idx = PYB11readwrite(doc="index of a 'representative' block")
    # alt_nodenum_vars = PYB11property(getterraw = """[](DBmultimesh& self) { std::vector<std::string> result;
    #                                                                  for (auto i = 0; i < self.nblocks; ++i) result.push_back(std::string(self.alt_nodenum_vars[i]));
    #                                                                  return result; }""")
    # alt_zonenum_vars = PYB11property(getterraw = """[](DBmultimesh& self) { std::vector<std::string> result;
    #                                                                  for (auto i = 0; i < self.nblocks; ++i) result.push_back(std::string(self.alt_zonenum_vars[i]));
    #                                                                  return result; }""")
    meshnames_alloc = PYB11readwrite(doc="original alloc of meshnames as string list")

#-------------------------------------------------------------------------------
class DBzonelist:
    "A silo zone list"
    ndims = PYB11readwrite(doc="Number of dimensions (2,3)")
    nzones = PYB11readwrite(doc="Number of zones in list")
    nshapes = PYB11readwrite(doc="Number of zone shapes")
    shapecnt = PYB11property(getterraw="[](DBzonelist& self) { return copy2py(self.shapecnt, self.nshapes); }",
                             setterraw="[](DBzonelist& self, py::list vals) { copy2c(self.shapecnt, vals); }",
                             doc="[nshapes] occurences of each shape")
    shapesize = PYB11property(getterraw="[](DBzonelist& self) { return copy2py(self.shapesize, self.nshapes); }",
                              setterraw="[](DBzonelist& self, py::list vals) { copy2c(self.shapesize, vals); }",
                              doc="[nshapes] Number of nodes per shape")
    shapetype = PYB11property(getterraw="[](DBzonelist& self) { return copy2py(self.shapetype, self.nshapes); }",
                              setterraw="[](DBzonelist& self, py::list vals) { copy2c(self.shapetype, vals); }",
                              doc="[nshapes] Type of shape")
    nodelist = PYB11property(getterraw="""[](DBzonelist& self) { size_t n = 0; 
                                                                 for (auto i = 0; i < self.nshapes; ++i) n += self.shapecnt[i]*self.shapesize[i];
                                                                 return copy2py(self.nodelist, n);
                                                               }""",
                             setterraw="[](DBzonelist& self, py::list vals) { copy2c(self.nodelist, vals); }",
                             doc="Sequent lst of nodes which comprise zones", returnpolicy="reference")
    lnodelist = PYB11readwrite(doc="Number of nodes in nodelist")
    origin = PYB11readwrite(doc="0' or '1'")
    min_index = PYB11readwrite(doc="Index of first real zone")
    max_index = PYB11readwrite(doc="Index of last real zone")
    zoneno = PYB11readwrite(doc="[nzones] zone number of each zone")
    gzoneno = PYB11readwrite(doc="[nzones] global zone number of each zone")
    gnznodtype = PYB11readwrite(doc="datatype for global node/zone ids")
    ghost_zone_labels = PYB11readwrite()

#-------------------------------------------------------------------------------
class DBphzonelist:
    "A silo polyhedral zone list"
    nfaces = PYB11readwrite(doc='Number of faces in facelist (aka "facetable")')
    nodecnt = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.nodecnt, self.nfaces); }",
                            setterraw="[](DBphzonelist& self, py::list vals) { copy2c(self.nodecnt, vals); }",
                            doc="Count of nodes in each face")
    lnodelist = PYB11readwrite(doc="Length of nodelist used to construct faces")
    nodelist = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.nodelist, self.lnodelist); }",
                             setterraw="[](DBphzonelist& self, py::list vals) { copy2c(self.nodelist, vals); }",
                             doc="List of nodes used in all faces")
    extface = PYB11readwrite(doc="boolean flag indicating if a face is external")
    nzones = PYB11readwrite(doc="Number of zones in this zonelist")
    facecnt = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.facecnt, self.nzones); }",
                            setterraw="[](DBphzonelist& self, py::list vals) { copy2c(self.facecnt, vals); }",
                            doc="Count of faces in each zone")
    lfacelist = PYB11readwrite(doc="Length of facelist used to construct zones")
    facelist = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.facelist, self.lfacelist); }",
                             setterraw="[](DBphzonelist& self, py::list vals) { copy2c(self.facelist, vals); }",
                             doc="List of faces used in all zones")
    origin = PYB11readwrite(doc="0' or '1'")
    lo_offset = PYB11readwrite(doc="Index of first non-ghost zone")
    hi_offset = PYB11readwrite(doc="Index of last non-ghost zone")
    zoneno = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.zoneno, self.nzones); }",
                           setterraw="[](DBphzonelist& self, py::list vals) { copy2c(self.zoneno, vals); }",
                           doc="[nzones] zone number of each zone")
    # gzoneno = PYB11property(getterraw="[](DBphzonelist& self) { return copy2py(self.gzoneno, self.nzones); }",
    #                         doc="[nzones] global zone number of each zone")
    gnznodtype = PYB11readwrite(doc="datatype for global node/zone ids")
    ghost_zone_labels = PYB11readwrite()

#-------------------------------------------------------------------------------
class DBfacelist:
    "A silo face list"
    ndims = PYB11readwrite(doc="Number of dimensions (2,3)")
    nfaces = PYB11readwrite(doc="Number of faces in list")
    origin = PYB11readwrite(doc="0' or '1'")
    nodelist = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.nodelist, self.lnodelist); }",
                             setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.nodelist, vals); }",
                             doc="Sequent list of nodes comprise faces")
    lnodelist = PYB11readwrite(doc="Number of nodes in nodelist")
    nshapes = PYB11readwrite(doc="Number of face shapes")
    shapecnt = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.shapecnt, self.nshapes); }",
                             setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.shapecnt, vals); }",
                             doc="[nshapes] Num of occurences of each shape")
    shapesize = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.shapesize, self.nshapes); }",
                             setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.shapesize, vals); }",
                              doc="[nshapes] Number of nodes per shape")
    ntypes = PYB11readwrite(doc="Number of face types")
    typelist = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.typelist, self.ntypes); }",
                             setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.typelist, vals); }",
                             doc="[ntypes] Type ID for each type")
    types = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.types, self.nfaces); }",
                          setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.types, vals); }",
                          doc="[nfaces] Type info for each face")
    nodeno = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.nodeno, self.lnodelist); }",
                           setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.nodeno, vals); }",
                           doc="[lnodelist] node number of each node")
    zoneno = PYB11property(getterraw="[](DBfacelist& self) { return copy2py(self.zoneno, self.nfaces); }",
                           setterraw="[](DBfacelist& self, py::list vals) { copy2c(self.zoneno, vals); }",
                           doc="[nfaces] Zone number for each face")

#-------------------------------------------------------------------------------
class DBedgelist:
    ndims = PYB11readwrite(doc="Number of dimensions (2,3)")
    nedges = PYB11readwrite(doc="Number of edges")
    edge_beg = PYB11property(getterraw="[](DBedgelist& self) { return copy2py(self.edge_beg, self.nedges); }",
                             setterraw="[](DBedgelist& self, py::list vals) { copy2c(self.edge_beg, vals); }")
    edge_end = PYB11property(getterraw="[](DBedgelist& self) { return copy2py(self.edge_end, self.nedges); }",
                             setterraw="[](DBedgelist& self, py::list vals) { copy2c(self.edge_end, vals); }")
    origin = PYB11readwrite(doc="0' or '1'")

#-------------------------------------------------------------------------------
class DBmultivar:
    id = PYB11readwrite(doc="Identifier for this object ")
    nvars = PYB11readwrite(doc="Number of variables  ")
    ngroups = PYB11readwrite(doc="Number of block groups in mesh")
    varnames = PYB11property(getterraw="[](DBmultivar& self) { return copy2py(self.varnames, self.nvars); }",
                             setterraw="[](DBmultivar& self, py::list vals) { copy2c(self.varnames, vals); }",
                             doc="Variable names")
    vartypes = PYB11property(getterraw="[](DBmultivar& self) { return copy2py(self.vartypes, self.nvars); }",
                             setterraw="[](DBmultivar& self, py::list vals) { copy2c(self.vartypes, vals); }",
                             doc="variable types")
    blockorigin = PYB11readwrite(doc="Origin (0 or 1) of block numbers")
    grouporigin = PYB11readwrite(doc="Origin (0 or 1) of group numbers")
    extentssize = PYB11readwrite(doc="size of each extent tuple")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    region_pnames = PYB11property(getterraw="[](DBmultivar& self) { return copy2py(self.region_pnames, self.ngroups); }",
                                  setterraw="[](DBmultivar& self, py::list vals) { copy2c(self.region_pnames, vals); }")
    mmesh_name = PYB11readwrite()
    tensor_rank = PYB11readwrite(doc="DB_VARTYPE_XXX")
    conserved = PYB11readwrite(doc="indicates if the variable should be conserved under various operations such as interp.")
    extensive = PYB11readwrite(doc="indicates if the variable reprsents an extensiv physical property (as opposed to intensive)")
    file_ns = PYB11readwrite(doc="namescheme for files (in lieu of meshnames)")
    block_ns = PYB11readwrite(doc="namescheme for block objects (in lieu of meshnames)")
    block_type = PYB11readwrite(doc="constant block type for all blocks (in lieu of meshtypes)")
    empty_list = PYB11property(getterraw="[](DBmultivar& self) { return copy2py(self.empty_list, self.empty_cnt); }",
                               setterraw="[](DBmultivar& self, py::list vals) { copy2c(self.empty_list, vals); }",
                               doc="list of empty block #'s (option for namescheme)")
    empty_cnt = PYB11readwrite(doc="size of empty list")
    repr_block_idx = PYB11readwrite(doc="index of a 'representative' block")
    missing_value = PYB11readwrite(doc="Value to indicate var data is invalid/missing")
    varnames_alloc = PYB11readwrite(doc="original alloc of varnames as string list")

#-------------------------------------------------------------------------------
class DBmultimat:
    id = PYB11readwrite(doc="Identifier for this object ")
    nmats = PYB11readwrite(doc="Number of materials  ")
    ngroups = PYB11readwrite(doc="Number of block groups in mesh")
    matnames = PYB11property(getterraw="[](DBmultimat& self) { return copy2py(self.matnames, self.nmats); }",
                             setterraw="[](DBmultimat& self, py::list vals) { copy2c(self.matnames, vals); }",
                             doc="names of constiuent DBmaterial objects")
    blockorigin = PYB11readwrite(doc="Origin (0 or 1) of block numbers")
    grouporigin = PYB11readwrite(doc="Origin (0 or 1) of group numbers")
    mixlens = PYB11property(getterraw="[](DBmultimat& self) { return copy2py(self.mixlens, self.nmats); }",
                            setterraw="[](DBmultimat& self, py::list vals) { copy2c(self.mixlens, vals); }",
                            doc="array of mixlen values in each mat")
    matcounts = PYB11property(getterraw="[](DBmultimat& self) { return copy2py(self.matcounts, self.ngroups); }",
                              setterraw="[](DBmultimat& self, py::list vals) { copy2c(self.matcounts, vals); }",
                              doc="counts of unique materials in each block")
    matlists = PYB11property(getterraw="""[](DBmultimat& self) { size_t nvals = 0;
                                                                 if (self.matcounts != NULL) {
                                                                   for (auto i = 0; i < self.ngroups; ++i) nvals += self.matcounts[i];
                                                                 }
                                                                 return copy2py(self.matlists, nvals);
                                                               }""",
                             setterraw="[](DBmultimat& self, py::list vals) { copy2c(self.matlists, vals); }",
                             doc="list of materials in each block")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    nmatnos = PYB11readwrite(doc="global number of materials over all pieces")
    matnos = PYB11property(getterraw="[](DBmultimat& self) { return copy2py(self.matnos, self.nmatnos); }",
                           setterraw="[](DBmultimat& self, py::list vals) { copy2c(self.matnos, vals); }",
                           doc="global list of material numbers")
    matcolors = PYB11property(getterraw = "[](DBmultimat& self) { return copy2py(self.matcolors, self.nmatnos); }",
                              setterraw = "[](DBmultimat& self, py::list vals) { copy2c(self.matcolors, vals); }",
                              doc = "optional colors for materials")
    material_names = PYB11property(getterraw = "[](DBmultimat& self) { return copy2py(self.material_names, self.nmatnos); }",
                                   setterraw = "[](DBmultimat& self, py::list vals) { copy2c(self.material_names, vals); }",
                                   doc = "optional names of the materials")
    allowmat0 = PYB11readwrite(doc='Flag to allow material "0"')
    mmesh_name = PYB11readwrite()
    file_ns = PYB11readwrite(doc="namescheme for files (in lieu of meshnames)")
    block_ns = PYB11readwrite(doc="namescheme for block objects (in lieu of meshnames)")
    empty_list = PYB11readwrite(doc="list of empty block #'s (option for namescheme)")
    empty_cnt = PYB11readwrite(doc="size of empty list")
    repr_block_idx = PYB11readwrite(doc="index of a 'representative' block")

#-------------------------------------------------------------------------------
class DBmultimatspecies:
    id = PYB11readwrite(doc="Identifier for this object ")
    nspec = PYB11readwrite(doc="Number of species  ")
    ngroups = PYB11readwrite(doc="Number of block groups in mesh")
    specnames = PYB11property(getterraw="[](DBmultimatspecies& self) { return copy2py(self.specnames, self.nspec); }",
                              setterraw="[](DBmultimatspecies& self, py::list vals) { copy2c(self.specnames, vals); }",
                              doc="Species object names")
    blockorigin = PYB11readwrite(doc="Origin (0 or 1) of block numbers")
    grouporigin = PYB11readwrite(doc="Origin (0 or 1) of group numbers")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    nmat = PYB11readwrite(doc="equiv. to nmatnos of a DBmultimat")
    nmatspec = PYB11readwrite(doc="equiv. to matnos of a DBmultimat")
    species_names = PYB11property(getterraw="[](DBmultimatspecies& self) { return copy2py(self.species_names, self.nspec); }",
                                  setterraw="[](DBmultimatspecies& self, py::list vals) { copy2c(self.species_names, vals); }",
                                  doc="optional names of the species")
    speccolors = PYB11property(getterraw="[](DBmultimatspecies& self) { return copy2py(self.speccolors, self.nspec); }",
                               setterraw="[](DBmultimatspecies& self, py::list vals) { copy2c(self.speccolors, vals); }",
                               doc="optional colors for species")
    file_ns = PYB11readwrite(doc="namescheme for files (in lieu of meshnames)")
    block_ns = PYB11readwrite(doc="namescheme for block objects (in lieu of meshnames)")
    empty_list = PYB11readwrite(doc="list of empty block #'s (option for namescheme)")
    empty_cnt = PYB11readwrite(doc="size of empty list")
    repr_block_idx = PYB11readwrite(doc="index of a 'representative' block")
    specnames_alloc = PYB11readwrite(doc="original alloc of matnames as string list")

#-------------------------------------------------------------------------------
class DBquadmesh:
    id = PYB11readwrite(doc="Identifier for this object")
    block_no = PYB11readwrite(doc="Block number for this mesh")
    group_no = PYB11readwrite(doc="Block group number for this mesh")
    name = PYB11readwrite(doc="Name associated with mesh")
    cycle = PYB11readwrite(doc="Problem cycle number")
    coord_sys = PYB11readwrite(doc="Cartesian, cylindrical, spherical")
    major_order = PYB11readwrite(doc="1 indicates row-major for multi-d arrays")
    stride = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.stride, 3); }",
                           setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.stride, vals); }",
                           doc="Offsets to adjacent elements")
    coordtype = PYB11readwrite(doc="Coord array type: collinear, non-collinear")
    facetype = PYB11readwrite(doc="Zone face type: rect, curv")
    planar = PYB11readwrite(doc="Sentinel: zones represent area or volume?")
    coords = PYB11property(getterraw = "[](DBquadmesh& self) { return copy2py(self.coords, 3, self.nnodes, self.datatype); }",
                           setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.coords, vals); }",
                           doc = "Mesh node coordinates")
    datatype = PYB11readwrite(doc="Type of coordinate arrays (double,float)")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    labels = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.labels, 3); }",
                           setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.labels, vals); }",
                           doc="Label associated with each dimension")
    units = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.units, 3); }",
                          setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.units, vals); }",
                          doc="Units for variable, e.g, 'mm/ms'")
    ndims = PYB11readwrite(doc="Number of computational dimensions")
    nspace = PYB11readwrite(doc="Number of physical dimensions")
    nnodes = PYB11readwrite(doc="Total number of nodes")
    dims = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.dims, 3); }",
                         setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.dims, vals); }",
                         doc="Number of nodes per dimension")
    origin = PYB11readwrite(doc="0' or '1'")
    min_index = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.min_index, 3); }",
                              setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.min_index, vals); }",
                              doc="Index in each dimension of 1st non-phoney")
    max_index = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.max_index, 3); }",
                              setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.max_index, vals); }",
                              doc="Index in each dimension of last non-phoney")
    base_index = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.base_index, 3); }",
                               setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.base_index, vals); }",
                               doc="Lowest real i,j,k value for this block")
    start_index = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.start_index, 3); }",
                                setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.start_index, vals); }",
                                doc="i,j,k values corresponding to original mesh")
    size_index = PYB11property(getterraw="[](DBquadmesh& self) { return copy2py(self.size_index, 3); }",
                               setterraw="[](DBquadmesh& self, py::list vals) { copy2c(self.size_index, vals); }",
                               doc="Number of nodes per dimension for original mesh")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    mrgtree_name = PYB11readwrite(doc="optional name of assoc. mrgtree object")
    ghost_node_labels = PYB11readwrite()
    ghost_zone_labels = PYB11readwrite()

#-------------------------------------------------------------------------------
class DBucdmesh:
    "A silo Unstructured Cell Data (UCD) Mesh"

    id = PYB11readwrite(doc="Identifier for this object")
    block_no = PYB11readwrite(doc="Block number for this mesh")
    group_no = PYB11readwrite(doc="Block group number for this mesh")
    name = PYB11readwrite(doc="Name associated with mesh")
    cycle = PYB11readwrite(doc="Problem cycle number")
    coord_sys = PYB11readwrite(doc="Coordinate system")
    topo_dim = PYB11readwrite(doc="Topological dimension. ")
    units = PYB11property(getterraw = "[](DBucdmesh& self) { return copy2py(self.units, 3); }",
                          setterraw="[](DBucdmesh& self, py::list vals) { copy2c(self.units, vals); }",
                          doc = "Units for variable, e.g, 'mm/ms'")
    labels = PYB11property(getterraw = "[](DBucdmesh& self) { return copy2py(self.labels, 3); }",
                           setterraw="[](DBucdmesh& self, py::list vals) { copy2c(self.labels, vals); }",
                           doc = "Label associated with each dimension")
    coords = PYB11property(getterraw = "[](DBucdmesh& self) { return copy2py(self.coords, 3, self.nnodes, self.datatype); }",
                           setterraw="[](DBucdmesh& self, py::list vals) { copy2c(self.coords, vals); }",
                           doc = "Mesh node coordinates")
    datatype = PYB11readwrite(doc="Type of coordinate arrays (double,float)")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    ndims = PYB11readwrite(doc="Number of computational dimensions")
    nnodes = PYB11readwrite(doc="Total number of nodes")
    origin = PYB11readwrite(doc="0' or '1'")
    faces = PYB11readwrite(doc="Data structure describing mesh faces", returnpolicy="reference")
    zones = PYB11readwrite(doc="Data structure describing mesh zones", returnpolicy="reference")
    edges = PYB11readwrite(doc="Data struct describing mesh edges", returnpolicy="reference")
    nodeno = PYB11readwrite(doc="nnodes] node number of each node")
    phzones = PYB11readwrite(doc="Data structure describing mesh zones")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    mrgtree_name = PYB11readwrite(doc="optional name of assoc. mrgtree object")
    tv_connectivity = PYB11readwrite()
    disjoint_mode = PYB11readwrite()
    gnznodtype = PYB11readwrite(doc="datatype for global node/zone ids")
    ghost_node_labels = PYB11readwrite()

#-------------------------------------------------------------------------------
# STL types
# vector_of_DBoptlist = PYB11_bind_vector("silo::DBoptlist_wrapper", opaque=True)

#-------------------------------------------------------------------------------
class DBquadvar:
    id = PYB11readwrite(doc="Identifier for this object")
    name = PYB11readwrite(doc="Name of variable")
    units = PYB11readwrite(doc="Units for variable, e.g, 'mm/ms'")
    label = PYB11readwrite(doc="Label (perhaps for editing purposes)")
    cycle = PYB11readwrite(doc="Problem cycle number")
    meshid = PYB11readwrite(doc="Identifier for associated mesh (Deprecated Sep2005)")
    vals = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.vals, self.nvals, self.nels, self.datatype); }",
                         setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.vals, vals, self.nvals, self.nels, self.datatype); }",
                         doc="Array of pointers to data arrays")
    datatype = PYB11readwrite(doc="Type of data pointed to by 'val'")
    nels = PYB11readwrite(doc="Number of elements in each array")
    nvals = PYB11readwrite(doc="Number of arrays pointed to by 'vals'")
    ndims = PYB11readwrite(doc="Rank of variable")
    dims = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.dims, 3); }",
                         setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.dims, vals); }",
                         doc = "Number of elements in each dimension")
    major_order = PYB11readwrite(doc="1 indicates row-major for multi-d arrays")
    stride = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.stride, 3); }",
                           setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.stride, vals); }",
                           doc = "Offsets to adjacent elements")
    min_index = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.min_index, 3); }",
                              setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.min_index, vals); }",
                              doc = "Index in each dimension of 1st non-phoney")
    max_index = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.max_index, 3); }",
                              setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.max_index, vals); }",
                              doc = "Index in each dimension of last non-phoney")
    origin = PYB11readwrite(doc="0' or '1'")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    mixvals = PYB11property(getterraw = "[](DBquadvar& self) { return copy2py(self.mixvals, self.nvals, self.mixlen, self.datatype); }",
                            setterraw = "[](DBquadvar& self, py::list vals) { copy2c(self.mixvals, vals); }",
                            doc = "nvals ptrs to data arrays for mixed zones")
    mixlen = PYB11readwrite(doc="Num of elmts in each mixed zone data array")
    use_specmf = PYB11readwrite(doc="Flag indicating whether to apply species mass fractions to the variable.")
    ascii_labels = PYB11readwrite(doc="Treat variable values as ASCII values by rounding to the nearest integer in the range [0, 255]")
    meshname = PYB11readwrite(doc="Name of associated mesh")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    conserved = PYB11readwrite(doc="indicates if the variable should be conserved under various operations such as interp.")
    extensive = PYB11readwrite(doc="indicates if the variable reprsents an extensive physical property (as opposed to intensive)")
    centering = PYB11readwrite(doc="explicit centering knowledge; should agree with alignment.")
    missing_value = PYB11readwrite(doc="Value to indicate var data is invalid/missing")

#-------------------------------------------------------------------------------
class DBucdvar:
    id = PYB11readwrite(doc="Identifier for this object")
    name = PYB11readwrite(doc="Name of variable")
    cycle = PYB11readwrite(doc="Problem cycle number")
    units = PYB11readwrite(doc="Units for variable, e.g, 'mm/ms'")
    label = PYB11readwrite(doc="Label (perhaps for editing purposes)")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    meshid = PYB11readwrite(doc="Identifier for associated mesh (Deprecated Sep2005)")
    vals = PYB11property(getterraw = "[](DBucdvar& self) { return copy2py(self.vals, self.nvals, self.nels, self.datatype); }",
                         setterraw = "[](DBucdvar& self, py::list vals) { copy2c(self.vals, vals); }",
                         doc = "Array of pointers to data arrays")
    datatype = PYB11readwrite(doc="Type of data pointed to by 'vals'")
    nels = PYB11readwrite(doc="Number of elements in each array")
    nvals = PYB11readwrite(doc="Number of arrays pointed to by 'vals'")
    ndims = PYB11readwrite(doc="Rank of variable")
    origin = PYB11readwrite(doc="0' or '1'")
    centering = PYB11readwrite(doc="Centering within mesh (nodal or zonal)")
    mixvals = PYB11property(getterraw = "[](DBucdvar& self) { return copy2py(self.mixvals, self.nvals, self.mixlen, self.datatype); }",
                            setterraw = "[](DBucdvar& self, py::list vals) { copy2c(self.mixvals, vals); }",
                           doc = "nvals ptrs to data arrays for mixed zones")
    mixlen = PYB11readwrite(doc="Num of elmts in each mixed zone data array")
    use_specmf = PYB11readwrite(doc="Flag indicating whether to apply species mass fractions to the variable.")
    ascii_labels = PYB11readwrite(doc="Treat variable values as ASCII values by rounding to the nearest integer in the range [0, 255]")
    meshname = PYB11readwrite(doc="Name of associated mesh")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    conserved = PYB11readwrite(doc="indicates if the variable should be conserved under various operations such as interp.")
    extensive = PYB11readwrite(doc="indicates if the variable reprsents an extensiv physical property (as opposed to intensive)")
    missing_value = PYB11readwrite(doc="Value to indicate var data is invalid/missing")

#-------------------------------------------------------------------------------
class DBmeshvar:
    "/*----------- Generic Mesh-Data Variable -----------*/"
    id = PYB11readwrite(doc="Identifier for this object")
    name = PYB11readwrite(doc="Name of variable")
    units = PYB11readwrite(doc="Units for variable, e.g, 'mm/ms'")
    label = PYB11readwrite(doc="Label (perhaps for editing purposes)")
    cycle = PYB11readwrite(doc="Problem cycle number")
    meshid = PYB11readwrite(doc="Identifier for associated mesh (Deprecated Sep2005)")
    vals = PYB11property(getterraw = "[](DBmeshvar& self) { return copy2py(self.vals, self.nvals, self.nels, self.datatype); }",
                         setterraw = "[](DBmeshvar& self, py::list vals) { copy2c(self.vals, vals); }",
                         doc = "Array of pointers to data arrays")
    datatype = PYB11readwrite(doc="Type of data pointed to by 'val'")
    nels = PYB11readwrite(doc="Number of elements in each array")
    nvals = PYB11readwrite(doc="Number of arrays pointed to by 'vals'")
    nspace = PYB11readwrite(doc="Spatial rank of variable")
    ndims = PYB11readwrite(doc="Rank of 'vals' array(s) (computatnl rank)")
    origin = PYB11readwrite(doc="0' or '1'")
    centering = PYB11readwrite(doc="Centering within mesh (nodal,zonal,other)")
    time = PYB11readwrite(doc="Problem time")
    dtime = PYB11readwrite(doc="Problem time, double data type")
    dims = PYB11property(getterraw = "[](DBmeshvar& self) { return copy2py(self.dims, 3); }",
                         setterraw = "[](DBmeshvar& self, py::list vals) { copy2c(self.dims, vals); }",
                         doc = "Number of elements in each dimension")
    major_order = PYB11readwrite(doc="1 indicates row-major for multi-d arrays")
    stride = PYB11property(getterraw = "[](DBmeshvar& self) { return copy2py(self.stride, 3); }",
                           setterraw = "[](DBmeshvar& self, py::list vals) { copy2c(self.stride, vals); }",
                           doc = "Offsets to adjacent elements")
    min_index = PYB11property(getterraw = "[](DBmeshvar& self) { return copy2py(self.min_index, 3); }",
                              setterraw = "[](DBmeshvar& self, py::list vals) { copy2c(self.min_index, vals); }",
                              doc = "Index in each dimension of 1st non-phoney")
    max_index = PYB11property(getterraw = "[](DBmeshvar& self) { return copy2py(self.max_index, 3); }",
                              setterraw = "[](DBmeshvar& self, py::list vals) { copy2c(self.max_index, vals); }",
                              doc = "Index in each dimension of last non-phoney")
    ascii_labels = PYB11readwrite(doc="Treat variable values as ASCII values by rounding to the nearest integer in the range [0, 255]")
    meshname = PYB11readwrite(doc="Name of associated mesh")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    conserved = PYB11readwrite(doc="indicates if the variable should be conserved under various operations such as interp.")
    extensive = PYB11readwrite(doc="indicates if the variable reprsents an extensiv physical property (as opposed to intensive)")
    missing_value = PYB11readwrite(doc="Value to indicate var data is invalid/missing")

#-------------------------------------------------------------------------------
class DBmaterial:
    "/*----------- Material Information -----------*/"
    id = PYB11readwrite(doc="Identifier")
    name = PYB11readwrite(doc="Name of this material information block")
    ndims = PYB11readwrite(doc="Rank of 'matlist' variable")
    origin = PYB11readwrite(doc="0' or '1'")
    dims = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(self.dims, 3); }",
                         setterraw = "[](DBmaterial& self, py::list vals) { copy2c(self.dims, vals); }",
                         doc = "Number of elements in each dimension")
    major_order = PYB11readwrite(doc="1 indicates row-major for multi-d arrays")
    stride = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(self.stride, 3); }",
                           setterraw = "[](DBmaterial& self, py::list vals) { copy2c(self.stride, vals); }",
                           doc = "Offsets to adjacent elements")
    nmat = PYB11readwrite(doc="Number of materials")
    matnos = PYB11property(getterraw="[](DBmaterial& self) { return copy2py(self.matnos, self.nmat); }",
                           setterraw="[](DBmaterial& self, py::list vals) { copy2c(self.matnos, vals); }",
                           doc="Array [nmat] of valid material numbers")
    matnames = PYB11property(getterraw="[](DBmaterial& self) { return copy2py(self.matnames, self.nmat); }",
                             setterraw="[](DBmaterial& self, py::list vals) { copy2c(self.matnames, vals); }",
                             doc="Array of material names")
    matlist = PYB11property(getterraw="""[](DBmaterial& self) { size_t nvals = 1;
                                                                for (auto i = 0; i < self.ndims; ++i) nvals *= self.dims[i];
                                                                return copy2py(self.matlist, nvals);
                                                              }""",
                            setterraw="[](DBmaterial& self, py::list vals) { copy2c(self.matlist, vals); }",
                            doc="Array[nzone] w/ mat. number or mix index")
    mixlen = PYB11readwrite(doc="Length of mixed data arrays (mix_xxx)")
    datatype = PYB11readwrite(doc="Type of volume-fractions (double,float)")
    mix_vf = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(static_cast<double*>(self.mix_vf), self.mixlen); }",
                           setterraw = "[](DBmaterial& self, py::list vals) { copy2c(static_cast<double*>(self.mix_vf), vals); }",
                           doc = "Array [mixlen] of volume fractions")
    mix_next = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(self.mix_next, self.mixlen); }",
                             setterraw = "[](DBmaterial& self, py::list vals) { copy2c(self.mix_next, vals); }",
                             doc = "Array [mixlen] of mixed data indices")
    mix_mat = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(self.mix_mat, self.mixlen); }",
                            setterraw = "[](DBmaterial& self, py::list vals) { copy2c(self.mix_mat, vals); }",
                            doc = "Array [mixlen] of mixed data indices")
    mix_zone = PYB11property(getterraw = "[](DBmaterial& self) { return copy2py(self.mix_zone, self.mixlen); }",
                             setterraw = "[](DBmaterial& self, py::list vals) { copy2c(self.mix_zone, vals); }",
                             doc = "Array [mixlen] of back pointers to mesh")
    meshname = PYB11readwrite(doc="Name of associated mesh")
    allowmat0 = PYB11readwrite(doc='Flag to allow material "0"')
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")

#-------------------------------------------------------------------------------
class DBmatspecies:
    "/*----------- Species Information -----------*/"
    id = PYB11readwrite(doc="Identifier")
    name = PYB11readwrite(doc="Name of this matspecies information block")
    matname = PYB11readwrite(doc="Name of material object with which material species object is associated.")
    nmat = PYB11readwrite(doc="Number of materials")
    nmatspec = PYB11property(getterraw = "[](DBmatspecies& self) { return copy2py(self.nmatspec, self.nmat); }",
                             setterraw = "[](DBmatspecies& self, py::list vals) { copy2c(self.nmatspec, vals); }",
                             doc = "Array of lngth nmat of the num of material species associated with each material.")
    ndims = PYB11readwrite(doc="Rank of 'speclist' variable")
    dims = PYB11property(getterraw="[](DBmatspecies& self) { return copy2py(self.dims, 3); }",
                         setterraw="[](DBmatspecies& self, py::list vals) { copy2c(self.dims, vals); }",
                         doc=" Number of elements in each dimension of the 'speclist' variable.")
    major_order = PYB11readwrite(doc="1 indicates row-major for multi-d arrays")
    stride = PYB11property(getterraw = "[](DBmatspecies& self) { return copy2py(self.stride, 3); }",
                           setterraw = "[](DBmatspecies& self, py::list vals) { copy2c(self.stride, vals); }",
                           doc = "Offsts to adjacent elmts in 'speclist'")
    nspecies_mf = PYB11readwrite(doc="Total number of species mass fractions.")
    species_mf = PYB11property(getterraw = "[](DBmatspecies& self) { return copy2py(static_cast<double*>(self.species_mf), self.nspecies_mf); }",
                               setterraw = "[](DBmatspecies& self, py::list vals) { copy2c(static_cast<double*>(self.species_mf), vals); }",
                             doc = "Array of length nspecies_mf of mass frations of the material species.")
    # int           *speclist;    /* Zone array of dimensions described by ndims
    #                              * and dims.  Each element of the array is an
    #                              * index into one of the species mass fraction
    #                              * arrays.  A positive value is the index in
    #                              * the species_mf array of the mass fractions
    #                              * of the clean zone's material species:
    #                              * species_mf[speclist[i]] is the mass fraction
    #                              * of the first species of material matlist[i]
    #                              * in zone i. A negative value means that the
    #                              * zone is a mixed zone and that the array
    #                              * mix_speclist contains the index to the
    #                              * species mas fractions: -speclist[i] is the
    #                              * index in the 'mix_speclist' array for zone
    #                              * i. */
    mixlen = PYB11readwrite(doc="Length of 'mix_speclist' array.")
    # int           *mix_speclist;  /* Array of lgth mixlen of 1-orig indices
    #                                * into the 'species_mf' array.
    #                                * species_mf[mix_speclist[j]] is the index
    #                                * in array species_mf' of the first of the
    #                                * mass fractions for material
    #                                * mix_mat[j]. */

    datatype = PYB11readwrite(doc="Datatype of mass fraction data.")
    guihide = PYB11readwrite(doc="Flag to hide from post-processor's GUI")
    # char         **specnames;   /* Array of species names; length is sum of nmatspec   */
    # char         **speccolors;  /* Array of species colors; length is sum of nmatspec */

#-------------------------------------------------------------------------------
# Module methods
@PYB11returnpolicy("reference")
@PYB11cppname("DBCreate_wrap")
def DBCreate(pathName = "std::string",
             mode = "int",
             target = "int",
             fileInfo = "std::string",
             fileType = "int"):
    "Create a silo database"
    return "DBfile*"

@PYB11returnpolicy("reference")
@PYB11cppname("DBOpen_wrap")
def DBOpen(name = "std::string",
           type = "int",
           mode = "int"):
    "Open a silo database"
    return "DBfile*"

@PYB11namespace("silo")
def DBClose():
    "Close a silo database"

@PYB11namespace("silo")
def DBMkDir():
    "Make a new path in a silo database"

@PYB11namespace("silo")
def DBSetDir():
    "Go to the given path in a silo database"

@PYB11namespace("silo")
def DBGetDir():
    "Get the current path in the silo database."

@PYB11namespace("silo")
def DBCpDir():
    "Copy a directory to another."

@PYB11namespace("silo")
def DBPutMultimesh():
    "Write a multimesh"

@PYB11namespace("silo")
@PYB11returnpolicy("take_ownership")
def DBGetMultimesh():
    "Read a multimesh"

@PYB11namespace("silo")
def DBPutMultimat():
    "Write a multimat"

def DBGetMultimat():
    "Get a multimat"

@PYB11namespace("silo")
def DBPutMultivar():
    "Write a multivar"

def DBGetMultivar():
    "Get a multivar"

@PYB11namespace("silo")
def DBPutMaterial():
    "Write a material"

def DBGetMaterial():
    "Get a material"

@PYB11namespace("silo")
def DBPutUcdmesh():
    "Write a UCD mesh"

@PYB11namespace("silo")
@PYB11returnpolicy("take_ownership")
def DBGetUcdmesh():
    "Read a UCD mesh"

@PYB11namespace("silo")
def DBPutQuadmesh():
    "Write a quad mesh"

def DBGetQuadmesh():
    "Get a quad mesh"

@PYB11namespace("silo")
def DBPutDefvars():
    "Write var definitions"

def DBGetDefvars():
    "Get var definitions"

def DBPutZonelist():
    "Write a zonelist"

def DBGetZonelist():
    "Get a zonelist"

@PYB11namespace("silo")
def DBPutZonelist2():
    "Write a zonelist"

@PYB11namespace("silo")
def DBPutPHZonelist():
    "Write a polyhedral zonelist"

def DBGetPHZonelist():
    "Get a polyhedral zonelist"

@PYB11namespace("silo")
def DBPutPointmesh():
    "Write a point mesh"

def DBGetPointmesh():
    "Get a point mesh"

@PYB11namespace("silo")
def DBAddRegion():
    "Add a region to a silo data base"

@PYB11namespace("silo")
def DBSetCwr():
    "set cwr?"

@PYB11namespace("silo")
def DBGetCwr():
    "Get the cwr"

@PYB11namespace("silo")
def DBPutMrgtree():
    "Write an mrg tree"

#-------------------------------------------------------------------------------
@PYB11template("T")
@PYB11namespace("silo")
def DBWrite():
    "Write a %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBWrite_vector():
    "Write a std::vector<%(T)s> to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBWrite_vector_of_vector():
    "Write a std::vector<std::vector<%(T)s>> to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutCompoundarray():
    "Write a compound array of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBReadVar():
    "Read a %(T)s from a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutUcdvar1():
    "Write a UCD scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutQuadvar1():
    "Write a quad mesh scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutPointvar1():
    "Write a point mesh scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutUcdvar():
    "Write a UCD mesh variable of %(T)s to a silo database."

def DBGetUcdvar():
    "Read a UCD mesh variable from a silo database."

def DBGetUcdvar():
    "Read a UCD mesh variable from a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutQuadvar():
    "Write a quad mesh variable of %(T)s to a silo database."

def DBGetQuadvar():
    "Read a quad mesh variable from a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutPointvar():
    "Write a point mesh variable of %(T)s to a silo database."

for d in ("int", "double"):
    exec('''
DBWrite_%(d)s = PYB11TemplateFunction(DBWrite, ("%(d)s",), pyname="DBWrite")
DBWrite_vector_%(d)s = PYB11TemplateFunction(DBWrite_vector, ("%(d)s",), pyname="DBWrite")
DBWrite_vector_of_vector_%(d)s = PYB11TemplateFunction(DBWrite_vector_of_vector, ("%(d)s",), pyname="DBWrite")
DBPutCompoundarray_%(d)s = PYB11TemplateFunction(DBPutCompoundarray, ("%(d)s", ), pyname="DBPutCompoundarray")
DBReadVar_%(d)s = PYB11TemplateFunction(DBReadVar, ("%(d)s",), pyname="DBReadVar")
DBPutUcdvar1_%(d)s = PYB11TemplateFunction(DBPutUcdvar1, ("%(d)s",), pyname="DBPutUcdvar1")
DBPutQuadvar1_%(d)s = PYB11TemplateFunction(DBPutQuadvar1, ("%(d)s",), pyname="DBPutQuadvar1")
DBPutPointvar1_%(d)s = PYB11TemplateFunction(DBPutPointvar1, ("%(d)s",), pyname="DBPutPointvar1")
''' % {"d" : d})

for d in ("Vector", "Tensor", "SymTensor"):
    for ndim in (2, 3):
        exec(('''
DBPutUcdvar_%(ndim)i = PYB11TemplateFunction(DBPutUcdvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutUcdvar")
DBPutQuadvar_%(ndim)i = PYB11TemplateFunction(DBPutQuadvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutQuadvar")
DBPutPointvar_%(ndim)i = PYB11TemplateFunction(DBPutPointvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutPointvar")
''') % {"ndim" : ndim,
        "d" : d})

#-------------------------------------------------------------------------------
# Taken from the silo.h file, expose the #define variables as module attributes.
DB_ZONETYPE_POLYHEDRON = PYB11attr()
DB_ZONETYPE_TET = PYB11attr()
DB_ZONETYPE_PYRAMID = PYB11attr()
DB_ZONETYPE_PRISM = PYB11attr()
DB_ZONETYPE_HEX = PYB11attr()

DB_NETCDF = PYB11attr()
DB_PDB = PYB11attr()
DB_TAURUS = PYB11attr()
DB_UNKNOWN = PYB11attr()
DB_DEBUG = PYB11attr()
DB_HDF5X = PYB11attr()
DB_PDBP = PYB11attr()

DB_HDF5_SEC2_OBSOLETE = PYB11attr()
DB_HDF5_STDIO_OBSOLETE = PYB11attr()
DB_HDF5_CORE_OBSOLETE = PYB11attr()
DB_HDF5_MPIO_OBSOLETE = PYB11attr()
DB_HDF5_MPIOP_OBSOLETE = PYB11attr()

DB_H5VFD_DEFAULT = PYB11attr()
DB_H5VFD_SEC2 = PYB11attr()
DB_H5VFD_STDIO = PYB11attr()
DB_H5VFD_CORE = PYB11attr()
DB_H5VFD_LOG = PYB11attr()
DB_H5VFD_SPLIT = PYB11attr()
DB_H5VFD_DIRECT = PYB11attr()
DB_H5VFD_FAMILY = PYB11attr()
DB_H5VFD_MPIO = PYB11attr()
DB_H5VFD_MPIP = PYB11attr()
DB_H5VFD_SILO = PYB11attr()

DB_FILE_OPTS_H5_DEFAULT_DEFAULT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SEC2 = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_STDIO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_CORE = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_LOG = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SPLIT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_DIRECT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_FAMILY = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_MPIO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_MPIP = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr()

DB_HDF5 = PYB11attr()
DB_HDF5_SEC2 = PYB11attr()
DB_HDF5_STDIO = PYB11attr()
DB_HDF5_CORE = PYB11attr()
DB_HDF5_LOG = PYB11attr()
DB_HDF5_SPLIT = PYB11attr()
DB_HDF5_DIRECT = PYB11attr()
DB_HDF5_FAMILY = PYB11attr()
DB_HDF5_MPIO = PYB11attr()
DB_HDF5_MPIOP = PYB11attr()
DB_HDF5_MPIP = PYB11attr()
DB_HDF5_SILO = PYB11attr()

DB_NFILES = PYB11attr()
DB_NFILTERS = PYB11attr()

DBAll = PYB11attr()
DBNone = PYB11attr()
DBCalc = PYB11attr()
DBMatMatnos = PYB11attr()
DBMatMatlist = PYB11attr()
DBMatMixList = PYB11attr()
DBCurveArrays = PYB11attr()
DBPMCoords = PYB11attr()
DBPVData = PYB11attr()
DBQMCoords = PYB11attr()
DBQVData = PYB11attr()
DBUMCoords = PYB11attr()
DBUMFacelist = PYB11attr()
DBUMZonelist = PYB11attr()
DBUVData = PYB11attr()
DBFacelistInfo = PYB11attr()
DBZonelistInfo = PYB11attr()
DBMatMatnames = PYB11attr()
DBUMGlobNodeNo = PYB11attr()
DBZonelistGlobZoneNo = PYB11attr()
DBMatMatcolors = PYB11attr()
DBCSGMBoundaryInfo = PYB11attr()
DBCSGMZonelist = PYB11attr()
DBCSGMBoundaryNames = PYB11attr()
DBCSGVData = PYB11attr()
DBCSGZonelistZoneNames = PYB11attr()
DBCSGZonelistRegNames = PYB11attr()
DBMMADJNodelists = PYB11attr()
DBMMADJZonelists = PYB11attr()
DBPMGlobNodeNo = PYB11attr()

DBObjectType = PYB11enum(("DB_INVALID_OBJECT",
                          "DB_QUADRECT",
                          "DB_QUADCURV",
                          "DB_QUADMESH",
                          "DB_QUADVAR",
                          "DB_UCDMESH",
                          "DB_UCDVAR",
                          "DB_MULTIMESH",
                          "DB_MULTIVAR",
                          "DB_MULTIMAT",
                          "DB_MULTIMATSPECIES",
                          "DB_MULTIBLOCKMESH",
                          "DB_MULTIBLOCKVAR",
                          "DB_MULTIMESHADJ",
                          "DB_MATERIAL",
                          "DB_MATSPECIES",
                          "DB_FACELIST",
                          "DB_ZONELIST",
                          "DB_EDGELIST",
                          "DB_PHZONELIST",
                          "DB_CSGZONELIST",
                          "DB_CSGMESH",
                          "DB_CSGVAR",
                          "DB_CURVE",
                          "DB_DEFVARS",
                          "DB_POINTMESH",
                          "DB_POINTVAR",
                          "DB_ARRAY",
                          "DB_DIR",
                          "DB_VARIABLE",
                          "DB_MRGTREE",
                          "DB_GROUPELMAP",
                          "DB_MRGVAR",
                          "DB_USERDEF"),
                         export_values = True)

DBdatatype = PYB11enum(("DB_INT",
                        "DB_SHORT",
                        "DB_LONG",
                        "DB_FLOAT",
                        "DB_DOUBLE",
                        "DB_CHAR",
                        "DB_LONG_LONG",
                        "DB_NOTYPE"),
                       export_values = True)

DB_CLOBBER = PYB11attr()
DB_NOCLOBBER = PYB11attr()

DB_READ = PYB11attr()
DB_APPEND = PYB11attr()

DB_LOCAL = PYB11attr()
DB_SUN3 = PYB11attr()
DB_SUN4 = PYB11attr()
DB_SGI = PYB11attr()
DB_RS6000 = PYB11attr()
DB_CRAY = PYB11attr()
DB_INTEL = PYB11attr()

DBOPT_FIRST = PYB11attr()
DBOPT_ALIGN = PYB11attr()
DBOPT_COORDSYS = PYB11attr()
DBOPT_CYCLE = PYB11attr()
DBOPT_FACETYPE = PYB11attr()
DBOPT_HI_OFFSET = PYB11attr()
DBOPT_LO_OFFSET = PYB11attr()
DBOPT_LABEL = PYB11attr()
DBOPT_XLABEL = PYB11attr()
DBOPT_YLABEL = PYB11attr()
DBOPT_ZLABEL = PYB11attr()
DBOPT_MAJORORDER = PYB11attr()
DBOPT_NSPACE = PYB11attr()
DBOPT_ORIGIN = PYB11attr()
DBOPT_PLANAR = PYB11attr()
DBOPT_TIME = PYB11attr()
DBOPT_UNITS = PYB11attr()
DBOPT_XUNITS = PYB11attr()
DBOPT_YUNITS = PYB11attr()
DBOPT_ZUNITS = PYB11attr()
DBOPT_DTIME = PYB11attr()
DBOPT_USESPECMF = PYB11attr()
DBOPT_XVARNAME = PYB11attr()
DBOPT_YVARNAME = PYB11attr()
DBOPT_ZVARNAME = PYB11attr()
DBOPT_ASCII_LABEL = PYB11attr()
DBOPT_MATNOS = PYB11attr()
DBOPT_NMATNOS = PYB11attr()
DBOPT_MATNAME = PYB11attr()
DBOPT_NMAT = PYB11attr()
DBOPT_NMATSPEC = PYB11attr()
DBOPT_BASEINDEX = PYB11attr()
DBOPT_ZONENUM = PYB11attr()
DBOPT_NODENUM = PYB11attr()
DBOPT_BLOCKORIGIN = PYB11attr()
DBOPT_GROUPNUM = PYB11attr()
DBOPT_GROUPORIGIN = PYB11attr()
DBOPT_NGROUPS = PYB11attr()
DBOPT_MATNAMES = PYB11attr()
DBOPT_EXTENTS_SIZE = PYB11attr()
DBOPT_EXTENTS = PYB11attr()
DBOPT_MATCOUNTS = PYB11attr()
DBOPT_MATLISTS = PYB11attr()
DBOPT_MIXLENS = PYB11attr()
DBOPT_ZONECOUNTS = PYB11attr()
DBOPT_HAS_EXTERNAL_ZONES = PYB11attr()
DBOPT_PHZONELIST = PYB11attr()
DBOPT_MATCOLORS = PYB11attr()
DBOPT_BNDNAMES = PYB11attr()
DBOPT_REGNAMES = PYB11attr()
DBOPT_ZONENAMES = PYB11attr()
DBOPT_HIDE_FROM_GUI = PYB11attr()
DBOPT_TOPO_DIM = PYB11attr()
DBOPT_REFERENCE = PYB11attr()
DBOPT_GROUPINGS_SIZE = PYB11attr()
DBOPT_GROUPINGS = PYB11attr()
DBOPT_GROUPINGNAMES = PYB11attr()
DBOPT_ALLOWMAT0 = PYB11attr()
DBOPT_MRGTREE_NAME = PYB11attr()
DBOPT_REGION_PNAMES = PYB11attr()
DBOPT_TENSOR_RANK = PYB11attr()
DBOPT_MMESH_NAME = PYB11attr()
DBOPT_TV_CONNECTIVITY = PYB11attr()
DBOPT_DISJOINT_MODE = PYB11attr()
DBOPT_MRGV_ONAMES = PYB11attr()
DBOPT_MRGV_RNAMES = PYB11attr()
DBOPT_SPECNAMES = PYB11attr()
DBOPT_SPECCOLORS = PYB11attr()
DBOPT_LLONGNZNUM = PYB11attr()
DBOPT_CONSERVED = PYB11attr()
DBOPT_EXTENSIVE = PYB11attr()
DBOPT_MB_FILE_NS = PYB11attr()
DBOPT_MB_BLOCK_NS = PYB11attr()
DBOPT_MB_BLOCK_TYPE = PYB11attr()
DBOPT_MB_EMPTY_LIST = PYB11attr()
DBOPT_MB_EMPTY_COUNT = PYB11attr()
DBOPT_LAST = PYB11attr()
  
DBOPT_H5_FIRST = PYB11attr()
DBOPT_H5_VFD = PYB11attr()
DBOPT_H5_RAW_FILE_OPTS = PYB11attr()
DBOPT_H5_RAW_EXTENSION = PYB11attr()
DBOPT_H5_META_FILE_OPTS = PYB11attr()
DBOPT_H5_META_EXTENSION = PYB11attr()
DBOPT_H5_CORE_ALLOC_INC = PYB11attr()
DBOPT_H5_CORE_NO_BACK_STORE = PYB11attr()
DBOPT_H5_META_BLOCK_SIZE = PYB11attr()
DBOPT_H5_SMALL_RAW_SIZE = PYB11attr()
DBOPT_H5_ALIGN_MIN = PYB11attr()
DBOPT_H5_ALIGN_VAL = PYB11attr()
DBOPT_H5_DIRECT_MEM_ALIGN = PYB11attr()
DBOPT_H5_DIRECT_BLOCK_SIZE = PYB11attr()
DBOPT_H5_DIRECT_BUF_SIZE = PYB11attr()
DBOPT_H5_LOG_NAME = PYB11attr()
DBOPT_H5_LOG_BUF_SIZE = PYB11attr()
DBOPT_H5_MPIO_COMM = PYB11attr()
DBOPT_H5_MPIO_INFO = PYB11attr()
DBOPT_H5_MPIP_NO_GPFS_HINTS = PYB11attr()
DBOPT_H5_SIEVE_BUF_SIZE = PYB11attr()
DBOPT_H5_CACHE_NELMTS = PYB11attr()
DBOPT_H5_CACHE_NBYTES = PYB11attr()
DBOPT_H5_CACHE_POLICY = PYB11attr()
DBOPT_H5_FAM_SIZE = PYB11attr()
DBOPT_H5_FAM_FILE_OPTS = PYB11attr()
DBOPT_H5_USER_DRIVER_ID = PYB11attr()
DBOPT_H5_USER_DRIVER_INFO = PYB11attr()
DBOPT_H5_SILO_BLOCK_SIZE = PYB11attr()
DBOPT_H5_SILO_BLOCK_COUNT = PYB11attr()
DBOPT_H5_SILO_LOG_STATS = PYB11attr()
DBOPT_H5_SILO_USE_DIRECT = PYB11attr()
DBOPT_H5_LAST = PYB11attr()

DB_TOP = PYB11attr()
DB_NONE = PYB11attr()
DB_ALL = PYB11attr()
DB_ABORT = PYB11attr()
DB_SUSPEND = PYB11attr()
DB_RESUME = PYB11attr()
DB_ALL_AND_DRVR = PYB11attr()

DB_ROWMAJOR = PYB11attr()
DB_COLMAJOR = PYB11attr()

DB_NOTCENT = PYB11attr()
DB_NODECENT = PYB11attr()
DB_ZONECENT = PYB11attr()
DB_FACECENT = PYB11attr()
DB_BNDCENT = PYB11attr()
DB_EDGECENT = PYB11attr()
DB_BLOCKCENT = PYB11attr()

DB_CARTESIAN = PYB11attr()
DB_CYLINDRICAL = PYB11attr()
DB_SPHERICAL = PYB11attr()
DB_NUMERICAL = PYB11attr()
DB_OTHER = PYB11attr()

DB_RECTILINEAR = PYB11attr()
DB_CURVILINEAR = PYB11attr()

DB_AREA = PYB11attr()
DB_VOLUME = PYB11attr()

DB_ON = PYB11attr()
DB_OFF = PYB11attr()

DB_ABUTTING = PYB11attr()
DB_FLOATING = PYB11attr()

DB_VARTYPE_SCALAR = PYB11attr()
DB_VARTYPE_VECTOR = PYB11attr()
DB_VARTYPE_TENSOR = PYB11attr()
DB_VARTYPE_SYMTENSOR = PYB11attr()
DB_VARTYPE_ARRAY = PYB11attr()
DB_VARTYPE_MATERIAL = PYB11attr()
DB_VARTYPE_SPECIES = PYB11attr()
DB_VARTYPE_LABEL = PYB11attr()

DBCSG_QUADRIC_G = PYB11attr()
DBCSG_SPHERE_PR = PYB11attr()
DBCSG_ELLIPSOID_PRRR = PYB11attr()
DBCSG_PLANE_G = PYB11attr()
DBCSG_PLANE_X = PYB11attr()
DBCSG_PLANE_Y = PYB11attr()
DBCSG_PLANE_Z = PYB11attr()
DBCSG_PLANE_PN = PYB11attr()
DBCSG_PLANE_PPP = PYB11attr()
DBCSG_CYLINDER_PNLR = PYB11attr()
DBCSG_CYLINDER_PPR = PYB11attr()
DBCSG_BOX_XYZXYZ = PYB11attr()
DBCSG_CONE_PNLA = PYB11attr()
DBCSG_CONE_PPA = PYB11attr()
DBCSG_POLYHEDRON_KF = PYB11attr()
DBCSG_HEX_6F = PYB11attr()
DBCSG_TET_4F = PYB11attr()
DBCSG_PYRAMID_5F = PYB11attr()
DBCSG_PRISM_5F = PYB11attr()

DBCSG_QUADRATIC_G = PYB11attr()
DBCSG_CIRCLE_PR = PYB11attr()
DBCSG_ELLIPSE_PRR = PYB11attr()
DBCSG_LINE_G = PYB11attr()
DBCSG_LINE_X = PYB11attr()
DBCSG_LINE_Y = PYB11attr()
DBCSG_LINE_PN = PYB11attr()
DBCSG_LINE_PP = PYB11attr()
DBCSG_BOX_XYXY = PYB11attr()
DBCSG_ANGLE_PNLA = PYB11attr()
DBCSG_ANGLE_PPA = PYB11attr()
DBCSG_POLYGON_KP = PYB11attr()
DBCSG_TRI_3P = PYB11attr()
DBCSG_QUAD_4P = PYB11attr()

DBCSG_INNER = PYB11attr()
DBCSG_OUTER = PYB11attr()
DBCSG_ON = PYB11attr()
DBCSG_UNION = PYB11attr()
DBCSG_INTERSECT = PYB11attr()
DBCSG_DIFF = PYB11attr()
DBCSG_COMPLIMENT = PYB11attr()
DBCSG_XFORM = PYB11attr()
DBCSG_SWEEP = PYB11attr()

DB_PREORDER = PYB11attr()
DB_POSTORDER = PYB11attr()
DB_FROMCWR = PYB11attr()

DB_ZONETYPE_BEAM = PYB11attr()

DB_ZONETYPE_POLYGON = PYB11attr()
DB_ZONETYPE_TRIANGLE = PYB11attr()
DB_ZONETYPE_QUAD = PYB11attr()
