#-------------------------------------------------------------------------------
# Helper to define methods for pickling/unpickling Spheral types.
# This can be important with pyMPI, since pyMPI relies on picking data types
# for communication!
#-------------------------------------------------------------------------------
import copyreg, pickle
from SpheralCompiledPackages import *
FacetedVolume1d = Box1d
FacetedVolume2d = Polygon
FacetedVolume3d = Polyhedron

#-------------------------------------------------------------------------------
# Vector
#-------------------------------------------------------------------------------
copyreg.pickle(Vector1d, lambda x: (Vector1d, (tuple(x))))
copyreg.pickle(Vector2d, lambda x: (Vector2d, (tuple(x))))
copyreg.pickle(Vector3d, lambda x: (Vector3d, (tuple(x))))

#-------------------------------------------------------------------------------
# Tensor
#-------------------------------------------------------------------------------
copyreg.pickle(Tensor1d, lambda x: (Tensor1d, (tuple(x))))
copyreg.pickle(Tensor2d, lambda x: (Tensor2d, (tuple(x))))
copyreg.pickle(Tensor3d, lambda x: (Tensor3d, (tuple(x))))

#-------------------------------------------------------------------------------
# SymTensor
#-------------------------------------------------------------------------------
def construct_SymTensor2d(xx, xy, yy, zz):
    return SymTensor2d(xx, xy, xy, yy, zz)

def construct_SymTensor3d(xx, xy, xz, yy, yz, zz):
    return SymTensor3d(xx, xy, xz,
                       xy, yy, yz,
                       xz, yz, zz)

copyreg.constructor(construct_SymTensor2d)
copyreg.constructor(construct_SymTensor3d)

copyreg.pickle(SymTensor1d, lambda x: (SymTensor1d, (tuple(x))))
copyreg.pickle(SymTensor2d, lambda x: (construct_SymTensor2d, (tuple(x))))
copyreg.pickle(SymTensor3d, lambda x: (construct_SymTensor3d, (tuple(x))))

#-------------------------------------------------------------------------------
# ThirdRankTensor
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor2d(*args):
    result = ThirdRankTensor2d()
    for i, x in enumerate(args):
        result[i] = x
    return result

def construct_ThirdRankTensor3d(*args):
    result = ThirdRankTensor3d()
    for i, x in enumerate(args):
        result[i] = x
    return result

copyreg.constructor(construct_ThirdRankTensor2d)
copyreg.constructor(construct_ThirdRankTensor3d)

copyreg.pickle(ThirdRankTensor1d, lambda x: (ThirdRankTensor1d, (tuple(x))))
copyreg.pickle(ThirdRankTensor2d, lambda x: (construct_ThirdRankTensor2d, (tuple(x))))
copyreg.pickle(ThirdRankTensor3d, lambda x: (construct_ThirdRankTensor3d, (tuple(x))))

#-------------------------------------------------------------------------------
# FactedVolume
#-------------------------------------------------------------------------------
def construct_Box1d(s):
    return unpackElementBox1d(s)

def construct_Polygon(s):
    return unpackElementPolygon(s)

def construct_Polyhedron(s):
    return unpackElementPolyhedron(s)

copyreg.constructor(construct_Box1d)
copyreg.constructor(construct_Polygon)
copyreg.constructor(construct_Polyhedron)

copyreg.pickle(Box1d, lambda x: (construct_Box1d, (packElementBox1d(x),)))
copyreg.pickle(Polygon, lambda x: (construct_Polygon, (packElementPolygon(x),)))
copyreg.pickle(Polyhedron, lambda x: (construct_Polyhedron, (packElementPolyhedron(x),)))

#------------------------------------------------------------------------------
# std::vectors of primitives
#------------------------------------------------------------------------------
def register_std_vector_primitive_for_pickle(T):
    vec = eval("vector_of_" + T)
    copyreg.constructor(vec)
    copyreg.pickle(vec, lambda x: (vec, (tuple(x),)))
for T in ("char", "unsigned", "ULL", "int", "float", "double", "string",
          "Vector1d", "Tensor1d", "SymTensor1d", "ThirdRankTensor1d", "FacetedVolume1d",
          "Vector2d", "Tensor2d", "SymTensor2d", "ThirdRankTensor2d", "FacetedVolume2d",
          "Vector3d", "Tensor3d", "SymTensor3d", "ThirdRankTensor3d", "FacetedVolume3d"):
    register_std_vector_primitive_for_pickle(T)

#------------------------------------------------------------------------------
# # std::vectors of general picklable types
# #------------------------------------------------------------------------------
# def register_std_vector_picklable_for_pickle(T):
#     vec = "vector_of_" + T
#     exec('''
# def make_std_vector_{T}(list_of_pickled_things):
#     return {vec}([pickle.loads(x) for x in list_of_pickled_things])
# def pickle_std_vector_{T}(xvec):
#     args = tuple([pickle.dumps(x) for x in xvec])
#     print("pickle_std_vector_{vec}: ", args)
#     return make_std_vector_{T}, (args,)
# '''.format(vec = vec,
#            T = T))
#     copyreg.constructor(eval("make_std_vector_" + T))
#     copyreg.pickle(eval(vec), eval("pickle_std_vector_" + T))
#     #copyreg.pickle(eval(vec), lambda xvec: (eval("make_std_vector_" + T), (tuple([pickle.dumps(x) for x in xvec]),)))
# for T in ("FacetedVolume1d",
#           "FacetedVolume2d",
#           "FacetedVolume3d"):
#     register_std_vector_picklable_for_pickle(T)
#     print("Registered vector_of_" + T)
