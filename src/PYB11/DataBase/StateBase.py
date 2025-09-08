#-------------------------------------------------------------------------------
# StateBase
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class StateBase:

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using FourthRankTensor = typename %(Dimension)s::FourthRankTensor;
    using FifthRankTensor = typename %(Dimension)s::FifthRankTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using KeyType = typename StateBase<%(Dimension)s>::KeyType;
    using FieldName = typename StateBase<%(Dimension)s>::FieldName;
    using MeshPtr = typename StateBase<%(Dimension)s>::MeshPtr;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def enroll(self, field="FieldBase<%(Dimension)s>&"):
        "Enroll a field to track"
        return "void"

    # @PYB11virtual
    # @PYB11pycppname("enroll")
    # def enroll_sharedptr(self, field="std::shared_ptr<FieldBase<%(Dimension)s>>&"):
    #     "Enroll a shared_ptr<Field> to track"
    #     return "void"

    @PYB11virtual
    @PYB11pycppname("enroll")
    def enroll_fieldlist(self, fieldList="FieldListBase<%(Dimension)s>&"):
        "Enroll a FieldList to track"
        return "void"

    @PYB11virtual
    def copyState(self):
        "Make an internal copy of all state referenced by this StateBase."
        return "void"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    #...........................................................................
    # Methods
    @PYB11const
    def registered(self, key="const KeyType&"):
        "Test if there is state registered corresponding to the Key"
        return "bool"

    @PYB11const
    @PYB11pycppname("registered")
    def registered1(self, field="const FieldBase<%(Dimension)s>&"):
        "Test if the specified Field is registered"
        return "bool"

    @PYB11const
    @PYB11pycppname("registered")
    def registered2(self, fieldList="const FieldListBase<%(Dimension)s>&"):
        "Test if the specified FieldList is registered"
        return "bool"

    @PYB11const
    def fieldNameRegistered(self, fieldName="const std::string&"):
        "Test if the a Field with the given name is registered"
        return "bool"

    def enrollMesh(self, meshPtr="MeshPtr"):
        "Enroll a mesh for tracking"
        return "void"

    # @PYB11pycppname("enroll")
    # def enroll_vec(self,
    #                key = "const std::string&",
    #                vec = "std::vector<Scalar>&"):
    #     "Enroll a vector<Scalar> using the given key"
    #     return "void"

    # @PYB11returnpolicy("reference_internal")
    # def array(self, key="const std::string&"):
    #     "Get the vector<double> associated with the given key"
    #     return "std::vector<Scalar>&"

    @PYB11const
    def keys(self):
        "The set of keys for state in the StateBase"
        return "std::vector<KeyType>"

    @PYB11const
    def fullFieldKeys(self):
        "The set of Field names (with NodeList mangling) for the state in the StateBase"
        return "std::vector<KeyType>"

    @PYB11const
    def fieldNames(self):
        "The set of unique Field names for the state in the StateBase (no NodeList mangling)"
        return "std::vector<KeyType>"

    @PYB11const
    def miscKeys(self):
        "The set of names for non-Fields in the StateBase"
        return "std::vector<KeyType>"

    def enrollConnectivityMap(self,
                              connectivityMapPtr = "std::shared_ptr<ConnectivityMap<%(Dimension)s>>"):
        "Enroll the ConnectivityMap"
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def connectivityMap(self):
        "Get the ConnectivityMap"
        return "const ConnectivityMap<%(Dimension)s>&"

    @PYB11const
    def meshRegistered(self):
        "Test if a mesh is registered"
        return "bool"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def mesh(self):
        "Get the currently registered Mesh"
        return "const Mesh<%(Dimension)s>&"

    def assign(self, rhs="const StateBase<%(Dimension)s>&"):
        "Set this StateBase's state equal to the other"
        return "void"

    @PYB11static
    def key(self, field="const FieldBase<%(Dimension)s>&"):
        "Construct a key for the given Field"
        return "KeyType"

    @PYB11static
    @PYB11pycppname("key")
    def key1(self, fieldList="const FieldListBase<%(Dimension)s>&"):
        "Construct a key for the given FieldList"
        return "KeyType"

    @PYB11static
    def buildFieldKey(self,
                      fieldName = "const KeyType&",
                      nodeListName = "const KeyType&"):
        "Hash a key based on the (fieldname, node list name)"
        return "KeyType"

    @PYB11static
    def splitFieldKey(self,
                      key = "const KeyType&",
                      fieldName = "KeyType&",
                      nodeListName = "KeyType&"):
        "Split the hashed key back into (field name, node list name) info"
        return "void"

    #...........................................................................
    # Template methods for getting Fields
    @PYB11template("Value")
    @PYB11returnpolicy("reference")
    @PYB11const
    def field(self,
              key = "const KeyType&",
              dummy = ("const %(Value)s&", "%(Value)s()")):
        "Return the %(Value)s field based on the key"
        return "Field<%(Dimension)s, %(Value)s>&"

    @PYB11template("Value")
    @PYB11const
    def fields(self,
               name = "const KeyType&",
               dummy = ("const %(Value)s&", "%(Value)s()"),
               allowNone = ("bool", "false")):
        "Return the %(Value)s FieldList based on the name"
        return "FieldList<%(Dimension)s, %(Value)s>"

    @PYB11template("Value")
    @PYB11const
    def allFields(self,
                  dummy = ("const %(Value)s&", "%(Value)s()")):
        "Return a set of all the %(Value)s Fields in the StateBase"
        return "std::vector<Field<%(Dimension)s, %(Value)s>*>"

    # unsignedField = PYB11TemplateMethod(field, "unsigned")
    # ULLField = PYB11TemplateMethod(field, "uint64_t")
    intField = PYB11TemplateMethod(field, "int")
    scalarField = PYB11TemplateMethod(field, "double")
    vectorField = PYB11TemplateMethod(field, "Vector")
    tensorField = PYB11TemplateMethod(field, "Tensor")
    symTensorField = PYB11TemplateMethod(field, "SymTensor")
    thirdRankTensorField = PYB11TemplateMethod(field, "ThirdRankTensor")
    fourthRankTensorField = PYB11TemplateMethod(field, "FourthRankTensor")
    fifthRankTensorField = PYB11TemplateMethod(field, "FifthRankTensor")
    facetedVolumeField = PYB11TemplateMethod(field, "FacetedVolume")
    vector_of_CellFaceFlagField = PYB11TemplateMethod(field, "std::vector<CellFaceFlag>")
    vector_of_doubleField = PYB11TemplateMethod(field, "std::vector<double>")
    RKCoefficientsField = PYB11TemplateMethod(field, "RKCoefficients<%(Dimension)s>")

    # unsignedFields = PYB11TemplateMethod(fields, "unsigned")
    # ULLFields = PYB11TemplateMethod(fields, "uint64_t")
    intFields = PYB11TemplateMethod(fields, "int")
    scalarFields = PYB11TemplateMethod(fields, "double")
    vectorFields = PYB11TemplateMethod(fields, "Vector")
    tensorFields = PYB11TemplateMethod(fields, "Tensor")
    symTensorFields = PYB11TemplateMethod(fields, "SymTensor")
    thirdRankTensorFields = PYB11TemplateMethod(fields, "ThirdRankTensor")
    fourthRankTensorFields = PYB11TemplateMethod(fields, "FourthRankTensor")
    fifthRankTensorFields = PYB11TemplateMethod(fields, "FifthRankTensor")
    facetedVolumeFields = PYB11TemplateMethod(fields, "FacetedVolume")
    vector_of_CellFaceFlagFields = PYB11TemplateMethod(fields, "std::vector<CellFaceFlag>")
    vector_of_doubleFields = PYB11TemplateMethod(fields, "std::vector<double>")
    RKCoefficientsFields = PYB11TemplateMethod(fields, "RKCoefficients<%(Dimension)s>")

    # allUnsignedFields = PYB11TemplateMethod(allFields, "unsigned")
    # allULLFields = PYB11TemplateMethod(allFields, "uint64_t")
    allIntFields = PYB11TemplateMethod(allFields, "int")
    allScalarFields = PYB11TemplateMethod(allFields, "double")
    allVectorFields = PYB11TemplateMethod(allFields, "Vector")
    allTensorFields = PYB11TemplateMethod(allFields, "Tensor")
    allSymTensorFields = PYB11TemplateMethod(allFields, "SymTensor")
    allThirdRankTensorFields = PYB11TemplateMethod(allFields, "ThirdRankTensor")
    allFourthRankTensorFields = PYB11TemplateMethod(allFields, "FourthRankTensor")
    allFifthRankTensorFields = PYB11TemplateMethod(allFields, "FifthRankTensor")
    allFacetedVolumeFields = PYB11TemplateMethod(allFields, "FacetedVolume")
    allVector_of_CellFaceFlagFields = PYB11TemplateMethod(allFields, "std::vector<CellFaceFlag>")
    allVector_of_doubleFields = PYB11TemplateMethod(allFields, "std::vector<double>")
    allRKCoefficientsFields = PYB11TemplateMethod(allFields, "RKCoefficients<%(Dimension)s>")

    #...........................................................................
    # enroll/get
    @PYB11template("Value")
    def enroll(self,
                  key = "const KeyType&",
                  thing = "%(Value)s&"):
        "Enroll a type of %(Value)s."
        return "void"

    @PYB11template("Value")
    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def get(self,
               key = "const KeyType&"):
        "Return a stored type of %(Value)s"
        return "%(Value)s&"

    enrollVectorVector = PYB11TemplateMethod(enroll, "std::vector<Vector>", pyname="enroll")
    getVectorVector = PYB11TemplateMethod(get, "std::vector<Vector>", pyname="get")

    #...........................................................................
    # assignFields
    @PYB11template("Value")
    def assignFields(self,
                     rhs = "const StateBase<%(Dimension)s>&",
                     name = "const std::string"):
        "Assign just the fields with the given name to those in another State object."
        return "void"

    assignFieldsScalar = PYB11TemplateMethod(assignFields, "double", pyname="assignFields")
    assignFieldsVector = PYB11TemplateMethod(assignFields, "Vector", pyname="assignFields")
    assignFieldsTensor = PYB11TemplateMethod(assignFields, "Tensor", pyname="assignFields")
    assignFieldsSymTensor = PYB11TemplateMethod(assignFields, "SymTensor", pyname="assignFields")
