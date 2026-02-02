#-------------------------------------------------------------------------------
# STL containers of primitives
#-------------------------------------------------------------------------------
from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

PYB11includes = ['"Geometry/Dimension.hh"',
                 '"RK/RKCoefficients.hh"',
                 "<vector>",
                 "<map>",
                 "<set>",
                 "<string>",
                 "<sstream>"]

# std::vector
vector_of_char     = PYB11_bind_vector("char", opaque=True, local=False)
vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True, local=False)
vector_of_ULL      = PYB11_bind_vector("uint64_t", opaque=True, local=False)
vector_of_int      = PYB11_bind_vector("int", opaque=True, local=False)
vector_of_float    = PYB11_bind_vector("float", opaque=True, local=False)
vector_of_double   = PYB11_bind_vector("double", opaque=True, local=False)
vector_of_string   = PYB11_bind_vector("std::string", opaque=True, local=False)

# std::vector<std::vector>
vector_of_vector_of_char     = PYB11_bind_vector("std::vector<char>", opaque=True, local=False)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True, local=False)
vector_of_vector_of_ULL      = PYB11_bind_vector("std::vector<uint64_t>", opaque=True, local=False)
vector_of_vector_of_int      = PYB11_bind_vector("std::vector<int>", opaque=True, local=False)
vector_of_vector_of_float    = PYB11_bind_vector("std::vector<float>", opaque=True, local=False)
vector_of_vector_of_double   = PYB11_bind_vector("std::vector<double>", opaque=True, local=False)
vector_of_vector_of_string   = PYB11_bind_vector("std::vector<std::string>", opaque=True, local=False)

# std::vector<pair<>>
vector_of_pair_double_double     = PYB11_bind_vector("std::pair<double, double>", opaque=True, local=False)
vector_of_pair_double_string     = PYB11_bind_vector("std::pair<double, std::string>", opaque=True, local=False)
vector_of_pair_unsigned_unsigned = PYB11_bind_vector("std::pair<unsigned, unsigned>", opaque=True, local=False)
vector_of_pair_ULL_ULL           = PYB11_bind_vector("std::pair<uint64_t, uint64_t>", opaque=True, local=False)
vector_of_pair_string_string     = PYB11_bind_vector("std::pair<std::string, std::string>", opaque=True, local=False)

# std::vector<tuple<>>
vector_of_tuple_int_int_int      = PYB11_bind_vector("std::tuple<int, int, int>", opaque=True, local=False)

# std::map
map_string_double = PYB11_bind_map("std::string", "double", opaque=True, local=False)
map_int_string    = PYB11_bind_map("int", "std::string", opaque=True, local=False)

@PYB11template("T1", "T2")
@PYB11namespace("std")
class pair:
    "A std::pair<%(T1)s, %(T2)s>"

    def pyinit(self):
        "Default constructor"
        return

    def pyinit(self,
               first = "%(T1)s",
               second = "%(T2)s"):
        "Construct with values (%(T1)s, %(T2)s)"
        return

    # @PYB11implementation('[](const std::pair<%(T1)s, %(T2)s>& self) { std::ostringstream os; os << "(" << self.first << ", " << self.second << ")"; return os.str(); }')
    # def __repr__(self):
    #     return

    first = PYB11readonly(doc="first value", returnpolicy="copy")
    second = PYB11readonly(doc="second value", returnpolicy="copy")

# std::pair
pair_double_string = PYB11TemplateClass(pair, template_parameters=("double", "std::string"))
pair_double_double = PYB11TemplateClass(pair, template_parameters=("double", "double"))

# RKCoefficients

for ndim in dims:
    Dimension = f"Spheral::Dim<{ndim}>"
    Vector = f"{Dimension}::Vector"
    exec(f'''
vector_of_RKCoefficients{ndim}d = PYB11_bind_vector("Spheral::RKCoefficients<{Dimension}>", opaque=True, local=False)
''')
