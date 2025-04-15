#ifndef __PBGWRAPS_SILOWRAPPERS__
#define __PBGWRAPS_SILOWRAPPERS__

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include "silo.h"

#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <memory>

namespace silo {

//------------------------------------------------------------------------------
template<typename T>
inline
py::list
copy2py(T* carray, const size_t nvals) {
  py::list result;
  for (auto i = 0; i < (int)nvals; ++i) result.append(carray[i]);
  return result;
}

inline
py::list
copy2py(char** carray, const size_t nvals) {
  py::list result;
  for (auto i = 0; i < (int)nvals; ++i) result.append(std::string(carray[i]));
  return result;
}

inline
py::list
copy2py(void** carray, const size_t nout, const size_t nin, const int dtype) {
  if (not (dtype == DB_FLOAT or dtype == DB_DOUBLE)) throw py::value_error("Require float or double type");
  py::list result;
  for (auto k = 0; k < (int)nout; ++k) {
    py::list row;
    if (dtype == DB_FLOAT) {
      auto* cvals = static_cast<float*>(carray[k]);
      for (auto i = 0; i < (int)nin; ++i) row.append(cvals[i]);
    } else {
      auto* cvals = static_cast<double*>(carray[k]);
      for (auto i = 0; i < (int)nin; ++i) row.append(cvals[i]);
    }
    result.append(row);
  }
  return result;
}

//------------------------------------------------------------------------------
template<typename T>
inline
void
copy2c(T* carray, py::list& vals) {
  if (carray != NULL) free(carray);
  const auto nvals = vals.size();
  carray = (T*) malloc(sizeof(T)*vals.size());
  for (auto i = 0; i < (int)nvals; ++i) carray[i] = vals[i].cast<T>();
}

inline
void
copy2c(char** carray, py::list& vals) {
  if (carray != NULL) delete [] carray;
  const auto nvals = vals.size();
  if (nvals > 0) {
    carray = new char*[nvals];
    for (auto i = 0; i < (int)nvals; ++i) {
      auto val = vals[i].cast<std::string>();
      carray[i] = new char[val.size() + 1];
      strcpy(carray[i], val.c_str());
    }
  }
}

inline
void
copy2c(void** carray, py::list& vals, const size_t nout, const size_t nin, const int dtype) {
  if (not (dtype == DB_FLOAT or dtype == DB_DOUBLE)) throw py::value_error("Require float or double type");
  if (vals.size() != nout) throw py::value_error("incorrectly sized list");
  for (auto k = 0; k < (int)nout; ++k) {
    auto rowvals = vals[k].cast<py::list>();
    if (rowvals.size() != nin) throw py::value_error("incorrectly sized inner list");
    if (carray[k] != NULL) free(carray[k]);
    if (dtype == DB_FLOAT) {
      carray[k] = malloc(sizeof(float)*nin);
      auto* cvals = static_cast<float*>(carray[k]);
      for (auto i = 0; i < (int)nin; ++i) cvals[i] = rowvals[i].cast<float>();
    } else {
      carray[k] = malloc(sizeof(double)*nin);
      auto* cvals = static_cast<double*>(carray[k]);
      for (auto i = 0; i < (int)nin; ++i) cvals[i] = rowvals[i].cast<double>();
    }
  }
}

//------------------------------------------------------------------------------
// Trait class to help mapping Spheral types to silo.
//------------------------------------------------------------------------------
template<typename T>
struct Spheral2Silo;

template<>
struct Spheral2Silo<double> {
  typedef double Value;
  static unsigned numElements() { return 1; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x;
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x = silovars[0][i];
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Vector> {
  typedef Spheral::Dim<2>::Vector Value;
  static unsigned numElements() { return 2; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.x(silovars[0][i]);
    x.y(silovars[1][i]);
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Tensor> {
  typedef Spheral::Dim<2>::Tensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.xx(silovars[0][i]); x.xy(silovars[1][i]);
    x.yx(silovars[2][i]); x.yy(silovars[3][i]);
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::SymTensor> {
  typedef Spheral::Dim<2>::SymTensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.xx(silovars[0][i]); x.xy(silovars[1][i]);
    x.yx(silovars[2][i]); x.yy(silovars[3][i]);
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Vector> {
  typedef Spheral::Dim<3>::Vector Value;
  static unsigned numElements() { return 3; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
    silovars[2][i] = x.z();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.x(silovars[0][i]);
    x.y(silovars[1][i]);
    x.z(silovars[2][i]);
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Tensor> {
  typedef Spheral::Dim<3>::Tensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.xx(silovars[0][i]); x.xy(silovars[1][i]); x.xz(silovars[2][i]);
    x.yx(silovars[3][i]); x.yy(silovars[4][i]); x.yz(silovars[5][i]);
    x.zx(silovars[6][i]); x.zy(silovars[7][i]); x.zz(silovars[8][i]);
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::SymTensor> {
  typedef Spheral::Dim<3>::SymTensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
  static void readElement(Value& x, double** silovars, const unsigned i) {
    x.xx(silovars[0][i]); x.xy(silovars[1][i]); x.xz(silovars[2][i]);
    x.yx(silovars[3][i]); x.yy(silovars[4][i]); x.yz(silovars[5][i]);
    x.zx(silovars[6][i]); x.zy(silovars[7][i]); x.zz(silovars[8][i]);
  }
};

//------------------------------------------------------------------------------
// Names
//------------------------------------------------------------------------------
// typedef ::DBfile DBfile;
// typedef ::DBoptlist DBoptlist;
// typedef ::DBmrgtree DBmrgtree;
// using ::DBMakeOptlist;
// using ::DBClearOptlist;
// using ::DBClearOption;
// using std::string;

//------------------------------------------------------------------------------
// Convert a std::string -> char*
//------------------------------------------------------------------------------
struct ConvertStringToCharStar {
  char* operator()(const std::string& x) {
    return const_cast<char*>(x.c_str());
  }
};

//------------------------------------------------------------------------------
// A trait class for for mapping types -> silo types.
//------------------------------------------------------------------------------
template<typename T> struct SiloTraits {};

template<>
struct SiloTraits<int> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_INT; }
};

template<>
struct SiloTraits<short> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_SHORT; }
};

template<>
struct SiloTraits<long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG; }
};

template<>
struct SiloTraits<float> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_FLOAT; }
};

template<>
struct SiloTraits<double> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_DOUBLE; }
};

template<>
struct SiloTraits<char> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_CHAR; }
};

template<>
struct SiloTraits<long long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG_LONG; }
};

//------------------------------------------------------------------------------
// Wrapper class to handle the memory management necessary with DBoptlist.
//------------------------------------------------------------------------------
struct DBoptlist_wrapper {
  DBoptlist* mOptlistPtr;
  std::vector<std::shared_ptr<void> > mCache;

  // Constructors.
  DBoptlist_wrapper(const int maxopts=1024):
    mOptlistPtr(DBMakeOptlist(maxopts)),
    mCache() {}

  // Destructor.
  ~DBoptlist_wrapper() {
    //ASSERT(DBFreeOptlist(mOptlistPtr) == 0);
  }

  // Generic functor definitions for adding and getting options.
  template<typename Value>
  struct AddOptionFunctor {
    int writeValue(DBoptlist_wrapper& optlist_wrapper,
                   const int option,
                   const Value& value) {
      std::shared_ptr<void> voidValue(new Value(value));
      optlist_wrapper.mCache.push_back(voidValue);
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, voidValue.get());
    }
    int writeVector(DBoptlist_wrapper& optlist_wrapper,
                    const int option,
                    const int option_size,
                    const std::vector<Value>& value) {
      DBoptlist_wrapper::AddOptionFunctor<int>().writeValue(optlist_wrapper, option_size, value.size());
      std::shared_ptr<void> voidValue(new std::vector<Value>(value));
      optlist_wrapper.mCache.push_back(voidValue);
      Value* frontPtr = &(((std::vector<Value>*) voidValue.get())->front());
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, frontPtr);
    }
  };

  template<typename Value>
  struct GetOptionFunctor {
    Value readValue(DBoptlist_wrapper& optlist_wrapper,
                    const int option) {
      return *((Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option));
    }
    std::vector<Value> readVector(DBoptlist_wrapper& optlist_wrapper,
                                  const int option,
                                  const int option_size) {
      const unsigned vecsize = DBoptlist_wrapper::GetOptionFunctor<int>().readValue(optlist_wrapper, option_size);
      Value* frontPtr = (Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
      return std::vector<Value>(frontPtr, frontPtr + vecsize);
    }
  };

  // Function definitions that use the functors.
  template<typename Value>
  int addOption(const int option,
                const Value& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeValue(*this, option, value);
  }

  template<typename Value>
  int addOption(const int option,
                const int option_size,
                const std::vector<Value>& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeVector(*this, option, option_size, value);
  }

  template<typename Value>
  Value getOption(const int option) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readValue(*this, option);
  }

  template<typename Value>
  std::vector<Value> getOption(const int option,
                               const int option_size) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readVector(*this, option, option_size);
  }
};
                              
//------------------------------------------------------------------------------
// Wrapper class to handle the memory managemnt necessary with DBmrgtree
//------------------------------------------------------------------------------
struct DBmrgtree_wrapper {
  DBmrgtree* mDBmrgtree;

  // Constructors.
  DBmrgtree_wrapper(int mesh_type,
                    int info_bits,
                    int max_children,
                    DBoptlist_wrapper optlist):
    mDBmrgtree(DBMakeMrgtree(mesh_type, info_bits, max_children, optlist.mOptlistPtr)) {}

  // Destructor.
  ~DBmrgtree_wrapper() {
    DBFreeMrgtree(mDBmrgtree);
  }

  // name
  std::string name() const { return (mDBmrgtree->name != NULL ? 
                                     std::string(mDBmrgtree->name) :
                                     std::string()); }
  void name(std::string val) { 
    mDBmrgtree->name = new char[val.length() + 1];
    strcpy(mDBmrgtree->name, val.c_str());
  }

  // src_mesh_name
  std::string src_mesh_name() const { return (mDBmrgtree->src_mesh_name != NULL ? 
                                              std::string(mDBmrgtree->src_mesh_name) :
                                              std::string()); }
  void src_mesh_name(std::string val) { 
    mDBmrgtree->src_mesh_name = new char[val.length() + 1];
    strcpy(mDBmrgtree->src_mesh_name, val.c_str());
  }

  // src_mesh_type
  int src_mesh_type() const { return mDBmrgtree->src_mesh_type; }
  void src_mesh_type(int val) { mDBmrgtree->src_mesh_type = val; }

  // type_info_bits
  int type_info_bits() const { return mDBmrgtree->type_info_bits; }
  void type_info_bits(int val) { mDBmrgtree->type_info_bits = val; }

  // num_nodes
  int num_nodes() const { return mDBmrgtree->num_nodes; }
  void num_nodes(int val) { mDBmrgtree->num_nodes = val; }
};
                              
//..............................................................................
// std::string specializations.
//..............................................................................
template<>
struct
DBoptlist_wrapper::AddOptionFunctor<std::string> {
  int
  writeValue(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const std::string& value) {
    std::shared_ptr<void> voidValue(new std::string(value));
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char*) ((std::string*) voidValue.get())->c_str());
  }
  int
  writeVector(DBoptlist_wrapper& optlist_wrapper,
              const int option,
              const int option_size,
              const std::vector<std::string>& value0) {
    VERIFY(optlist_wrapper.addOption<int>(option_size, value0.size()) == 0);
    std::shared_ptr<void> voidValue(new std::vector<std::string>(value0));
    auto value = static_cast<std::vector<std::string>*>(voidValue.get());
    std::shared_ptr<void> voidChar(new char*[value->size()]);
    char** charArray = (char**) voidChar.get();
    for (auto k = 0; k < (int)value->size(); ++k) {
      charArray[k] = new char[(*value)[k].size() + 1];
      strcpy(charArray[k], (*value)[k].c_str());
    }
    optlist_wrapper.mCache.push_back(voidChar);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, charArray);
  }
};

template<>
struct
DBoptlist_wrapper::GetOptionFunctor<std::string> {
  std::string
  readValue(DBoptlist_wrapper& optlist_wrapper,
            const int option) {
    char* result = (char*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    return std::string(result);
  }
  std::vector<std::string>
  readVector(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const int option_size) {
    const int resultsize = optlist_wrapper.getOption<int>(option_size);
    VERIFY(resultsize > 0);
    char** chararray = (char**) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    std::vector<std::string> result(chararray, chararray + resultsize);
    // for (unsigned k = 0;  k != resultsize; ++k) result.push_back(std::string(chararray[k]));
    VERIFY((int)result.size() == resultsize);
    return result;
  }
};

//------------------------------------------------------------------------------
// DBCreate
//------------------------------------------------------------------------------
inline
DBfile*
DBCreate_wrap(std::string pathName,
              int mode,
              int target,
              std::string fileInfo,
              int fileType) {
  auto* result  = DBCreate(pathName.c_str(), mode, target, fileInfo.c_str(), fileType);
  VERIFY2(result != nullptr, "Error creating file " << pathName);
  return result;
  // return DBCreate(pathName.c_str(), mode, target, fileInfo.c_str(), fileType);
}

//------------------------------------------------------------------------------
// DBOpen
//------------------------------------------------------------------------------
inline
DBfile*
DBOpen_wrap(std::string pathName,
              int type,
              int mode) {
  return DBOpen(pathName.c_str(), type, mode);
}

//------------------------------------------------------------------------------
// DBMakeMrgtree
//------------------------------------------------------------------------------
// inline
// DBmrgtree*
// DBMakeMrgtree_wrap(int mesh_type,
//                    int info_bits,
//                    int max_children,
//                    DBoptlist_wrapper& optlist) {
//   return DBMakeMrgtree(mesh_type,
//                        info_bits,
//                        max_children,
//                        optlist.mOptlistPtr);
// }

//------------------------------------------------------------------------------
// DBFreeMrgtree
//------------------------------------------------------------------------------
// inline
// void
// DBFreeMrgtree_wrap(DBmrgtree& tree) {
//   DBFreeMrgtree(&tree);
// }

//------------------------------------------------------------------------------
// DBClose
//------------------------------------------------------------------------------
inline
int
DBClose(DBfile& file) {
  return DBClose(&file);
}

//------------------------------------------------------------------------------
// DBMkDir
//------------------------------------------------------------------------------
inline
int
DBMkDir(DBfile& file,
        std::string dirname) {
  return DBMkDir(&file, dirname.c_str());
}

//------------------------------------------------------------------------------
// DBSetDir
//------------------------------------------------------------------------------
inline
int
DBSetDir(DBfile& file,
         std::string dirname) {
  return DBSetDir(&file, dirname.c_str());
}

//------------------------------------------------------------------------------
// DBGetDir
//------------------------------------------------------------------------------
inline
std::string
DBGetDir(DBfile& file) {
  char result[256];
  auto valid = DBGetDir(&file, result);
  VERIFY2(valid == 0, "Silo ERROR: unable to fetch directory name.");
  return std::string(result);
}

//------------------------------------------------------------------------------
// DBCpDir
//------------------------------------------------------------------------------
inline
int
DBCpDir(DBfile& srcFile,
        std::string srcDir,
        DBfile& dstFile,
        std::string dstDir) {
  return DBCpDir(&srcFile, srcDir.c_str(), &dstFile, dstDir.c_str());
}

//------------------------------------------------------------------------------
// DBWrite
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite(DBfile& file,
        std::string varname,
        T& var) {
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &var, 
                 &(SiloTraits<T>::dims()).front(),
                 SiloTraits<T>::dims().size(),
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBWrite py::sequence
//------------------------------------------------------------------------------
// template<typename T>
// inline
// int
// DBWrite_sequence(DBfile& file,
//                  std::string varname,
//                  py::sequence& var) {
//   const auto n = var.size();
//   auto dims = std::vector<int>(1, n);
//   std::vector<T> vals(n);
//   for (auto i = 0; i < n; ++i) vals[i] = var[i].cast<T>();
//   return DBWrite(&file,
//                  varname.c_str(),
//                  (void*) &vals.front(), 
//                  &dims.front(),
//                  1,
//                  SiloTraits<T>::datatype());
// }

//------------------------------------------------------------------------------
// DBWrite vector<T>
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite_vector(DBfile& file,
               std::string varname,
               std::vector<T>& var) {
  auto dims = std::vector<int>(1, var.size());
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &var.front(), 
                 &dims.front(),
                 1,
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBWrite vector<vector<T>>
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite_vector_of_vector(DBfile& file,
                         std::string varname,
                         std::vector<std::vector<T>>& var) {
  auto ndims = var.size();
  auto dims = std::vector<int>(ndims);
  std::vector<T> varlinear;
  for (auto i = 0; i < (int)ndims; ++i) {
    dims[i] = var[i].size();
    auto istart = varlinear.size();
    varlinear.resize(varlinear.size() + dims[i]);
    std::copy(var[i].begin(), var[i].end(), varlinear.begin() + istart);
  }
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &varlinear.front(),
                 &dims.front(),
                 ndims,
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBReadVar
//------------------------------------------------------------------------------
template<typename T>
inline
T
DBReadVar(DBfile& file,
          std::string varname) {
  T result;
  DBReadVar(&file, varname.c_str(), (void*) &result);
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultimesh
//------------------------------------------------------------------------------
inline
int
DBPutMultimesh(DBfile& file,
               std::string name,
               std::vector<std::string>& meshNames,
               std::vector<int>& meshTypes,
               DBoptlist_wrapper& optlist) {

  // Pre-conditions.
  VERIFY2(meshNames.size() == meshTypes.size(), "meshNames and meshTypes must be same length:  " << meshNames.size() << " != " << meshTypes.size());

  // Convert names to char*.
  std::vector<char*> meshNames1;
  std::transform(meshNames.begin(), meshNames.end(), std::back_inserter(meshNames1),
                 ConvertStringToCharStar());
  CHECK(meshNames1.size() == meshNames.size());

  // Do the deed.
  const int result = DBPutMultimesh(&file, 
                                    name.c_str(), 
                                    meshNames.size(),
                                    &meshNames1.front(), 
                                    &meshTypes.front(),
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBGetMultimesh
//------------------------------------------------------------------------------
inline
DBmultimesh*
DBGetMultimesh(DBfile& file,
               std::string name) {
  return DBGetMultimesh(&file, name.c_str());
}

//------------------------------------------------------------------------------
// DBPutMultimat
//------------------------------------------------------------------------------
inline
int
DBPutMultimat(DBfile& file,
              std::string name,
              std::vector<std::string>& matNames,
              DBoptlist_wrapper& optlist) {

  // Convert names to char*.
  std::vector<char*> matNames1;
  std::transform(matNames.begin(), matNames.end(), std::back_inserter(matNames1),
                 ConvertStringToCharStar());
  CHECK(matNames1.size() == matNames.size());

  // Do the deed.
  const int result = DBPutMultimat(&file, 
                                   name.c_str(), 
                                   matNames.size(),
                                   &matNames1.front(), 
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutCompoundarray
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutCompoundarray(DBfile& file,
                   std::string name,
                   std::vector<std::string>& elemNames,
                   std::vector<std::vector<T> >& values,
                   DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(elemNames.size() == values.size());

  // Convert names to char*.
  std::vector<char*> elemNames1;
  std::transform(elemNames.begin(), elemNames.end(), std::back_inserter(elemNames1),
                 ConvertStringToCharStar());
  CHECK(elemNames1.size() == elemNames.size());

  // Read the sizes of each array.
  std::vector<int> elemLengths;
  elemLengths.reserve(values.size());
  for (unsigned k = 0; k != values.size(); ++k) elemLengths.push_back(values[k].size());
  const unsigned numValues = std::accumulate(elemLengths.begin(), elemLengths.end(), 0);

  // Flatten the values to a single arrray.
  std::vector<T> flatValues;
  flatValues.reserve(numValues);
  for (unsigned k = 0; k != values.size(); ++k) {
    for (unsigned j = 0; j != values[k].size(); ++j) {
      flatValues.push_back(values[k][j]);
    }
  }
  CHECK(flatValues.size() == numValues);

  // Do the deed.
  CHECK(elemNames.size() == elemLengths.size());
  const int result = DBPutCompoundarray(&file, 
                                        name.c_str(), 
                                        &elemNames1.front(),
                                        &elemLengths.front(),
                                        elemNames1.size(),
                                        &flatValues.front(),
                                        numValues,
                                        SiloTraits<T>::datatype(),
                                        optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultivar
//------------------------------------------------------------------------------
inline
int
DBPutMultivar(DBfile& file,
              std::string name,
              std::vector<std::string>& varNames,
              std::vector<int>& varTypes,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned numVars = varNames.size();
  VERIFY(varTypes.size() == numVars);

  // Convert names to char*.
  std::vector<char*> varNames1;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(varNames1),
                 ConvertStringToCharStar());
  CHECK(varNames1.size() == varNames.size());

  // Do the deed.
  const int result = DBPutMultivar(&file, 
                                   const_cast<char*>(name.c_str()), 
                                   numVars,
                                   &varNames1.front(),
                                   &varTypes.front(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMaterial
//------------------------------------------------------------------------------
inline
int
DBPutMaterial(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<int>& matnos,
              std::vector<int>& matlist,
              std::vector<int>  dims,
              std::vector<int>& mix_next,
              std::vector<int>& mix_mat,
              std::vector<int>& mix_zone,
              std::vector<double>& mix_vf,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nmat = matnos.size();
  const unsigned numMix = mix_next.size();
  VERIFY(mix_mat.size() == numMix);
  VERIFY(mix_zone.size() == numMix);
  VERIFY(mix_vf.size() == numMix);

  // If dims is empty, set it as a 1D list based on matlist.
  if (dims.empty()) dims = std::vector<int>(1, matlist.size());

  // If mix info is empty, we need to be sure and point at valid empty arrays
  int *mix_next_begin, *mix_mat_begin, *mix_zone_begin;
  double *mix_vf_begin;
  if (numMix == 0) {
    mix_next_begin = nullptr;
    mix_mat_begin = nullptr;
    mix_zone_begin = nullptr;
    mix_vf_begin = nullptr;
  } else {
    mix_next_begin = &mix_next.front();
    mix_mat_begin = &mix_mat.front();
    mix_zone_begin = &mix_zone.front();
    mix_vf_begin = &mix_vf.front();
  }

  // Do the deed.
  const int result = DBPutMaterial(&file, 
                                   name.c_str(),
                                   meshName.c_str(),
                                   nmat,
                                   &matnos.front(),
                                   &matlist.front(),
                                   &dims.front(),
                                   dims.size(),
                                   mix_next_begin,
                                   mix_mat_begin,
                                   mix_zone_begin,
                                   mix_vf_begin,
                                   numMix,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutUcdmesh
//------------------------------------------------------------------------------
inline
int
DBPutUcdmesh(DBfile& file,
             std::string name,
             std::vector<std::vector<double> >& coords,
             int nzones,
             std::string zonel_name,
             std::string facel_name,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned nnodes = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == nnodes);
  VERIFY(nzones > 0);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[nnodes];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Convert strings to char*.
  char* zonel_name1 = (zonel_name == "NULL") ? NULL : const_cast<char*>(zonel_name.c_str());
  char* facel_name1 = (facel_name == "NULL") ? NULL : const_cast<char*>(facel_name.c_str());

  // Do the deed.
  const int result = DBPutUcdmesh(&file, 
                                  name.c_str(),
                                  ndims,
                                  NULL,
                                  coordPtrs,
                                  nnodes,
                                  nzones,
                                  zonel_name1,
                                  facel_name1,
                                  SiloTraits<double>::datatype(),
                                  optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBGetUcdmesh
//------------------------------------------------------------------------------
inline
DBucdmesh*
DBGetUcdmesh(DBfile& file,
             std::string name) {
  return DBGetUcdmesh(&file, name.c_str());
}

//------------------------------------------------------------------------------
// DBPutQuadmesh
// Note we assume just the unique (x,y,z) coordinates are provided, but we
// replicate them here for the silo file writing.
//------------------------------------------------------------------------------
inline
int
DBPutQuadmesh(DBfile& file,
              std::string name,
              std::vector<std::vector<double> >& coords,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const auto ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);

  // Number of nodes in each dimension.
  auto nxnodes = coords[0].size();
  auto nxynodes = nxnodes*coords[1].size();
  std::vector<int> meshdims(ndims);
  auto nnodes = 1;
  for (auto k = 0; k < (int)ndims; ++k) {
    meshdims[k] = coords[k].size();
    nnodes *= coords[k].size();
  }

  // We need the C-stylish pointers to the coordinates.
  // This is where we flesh out to the nnodes number of values too.
  double** coordPtrs = new double*[ndims];
  for (auto k = 0; k < (int)ndims; ++k) coordPtrs[k] = new double[nnodes];
  for (auto inode = 0; inode < nnodes; ++inode) {
    const size_t index[3] = {inode % nxnodes,
                             (inode % nxynodes) / nxnodes,
                             inode / nxynodes};
    for (auto k = 0; k < (int)ndims; ++k) coordPtrs[k][inode] = coords[k][index[k]];
  }

  // Do the deed.
  const int result = DBPutQuadmesh(&file,                            // dbfile
                                   name.c_str(),                     // name
                                   NULL,                             // coordnames
                                   coordPtrs,                        // coords
                                   &meshdims[0],                     // dims
                                   ndims,                            // ndims
                                   SiloTraits<double>::datatype(),   // datatype
                                   DB_NONCOLLINEAR,                  // coordtype
                                   optlist.mOptlistPtr);             // optlist

  // That's it.
  for (auto k = 0; k < (int)ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutDefvars
//------------------------------------------------------------------------------
inline
int
DBPutDefvars(DBfile& file,
             std::string name,
             std::vector<std::string>& varNames,
             std::vector<int>& varTypes,
             std::vector<std::string>& varDefs,
             std::vector<DBoptlist_wrapper*>& optlists) {

  // Preconditions.
  const unsigned ndefs = varNames.size();
  VERIFY(varNames.size() == ndefs and
         varTypes.size() == ndefs and
         varDefs.size() == ndefs and
         optlists.size() == ndefs);

  // Convert names to char*'s
  std::vector<char*> names, defns;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(names), ConvertStringToCharStar());
  std::transform(varDefs.begin(), varDefs.end(), std::back_inserter(defns), ConvertStringToCharStar());
  CHECK(names.size() == varNames.size());
  CHECK(defns.size() == varDefs.size());
  
  // Copy the optlists to an array of pointers.
  std::vector<DBoptlist*> optlistptrs;
  for (std::vector<DBoptlist_wrapper*>::iterator itr = optlists.begin();
       itr != optlists.end();
       ++itr) optlistptrs.push_back((*itr)->mOptlistPtr);
  CHECK(optlistptrs.size() == ndefs);

  return DBPutDefvars(&file,
                      name.c_str(),
                      ndefs,
                      &names.front(),
                      &varTypes.front(),
                      &defns.front(),
                      &optlistptrs.front());
}

//------------------------------------------------------------------------------
// DBPutUcdvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar(DBfile& file,
            std::string name,
            std::string meshName,
            std::vector<T>& values,
            std::vector<T>& mixValues,
            int centering,
            DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();
  int mixlen = mixValues.size();
  std::vector<char*> varnames;
  for (i = 0; (int)i != nvars; ++i) {
    varnames.push_back(const_cast<char*>((name + "_").c_str()));
    sprintf(varnames.back(), "%i", i);
  }

  // Build the sub-variables.
  double** vars = new double*[nvars];
  double** mixvars = new double*[nvars];
  for (i = 0; (int)i != nvars; ++i) {
    vars[i] = new double[nels];
    mixvars[i] = new double[mixlen];
  }
  for (j = 0; (int)j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);
  for (j = 0; (int)j != mixlen; ++j) Spheral2Silo<T>::copyElement(mixValues[j], mixvars, j);

  const int result = DBPutUcdvar(&file,
                                 name.c_str(),
                                 meshName.c_str(),
                                 nvars,
                                 &varnames.front(),
                                 (void*) vars,
                                 nels,
                                 (void*) mixvars,
                                 mixlen,
                                 SiloTraits<double>::datatype(),
                                 centering,
                                 optlist.mOptlistPtr);

  // That's it.
  for (i = 0; (int)i != nvars; ++i) {
    delete[] vars[i];
    delete[] mixvars[i];
  }
  delete[] vars;
  delete[] mixvars;
  return result;
}

// //------------------------------------------------------------------------------
// // DBReadUcdvar
// // We assume here that the underlying element type is double.
// //------------------------------------------------------------------------------
// template<typename T>
// inline
// void
// DBReadUcdvar(DBfile& file,
//              std::string name,
//              std::vector<T>& values,
//              std::vector<T>& mixValues) {

//   auto ucdvar = DBGetUcdvar(&file,
//                             name.c_str());
  
//   // Read into the resulting arrays.
//   const auto element_size = Spheral2Silo<T>::numElements();
//   const auto nvals = ucdvar->nels/element_size;
//   const auto nmixvals = ucdvar->mixlen/element_size;
//   values.resize(nvals);
//   mixValues.resize(nmixvals);
//   for (auto i = 0; i < nvals; ++i)    Spheral2Silo<T>::readElement(values[i], ucdvar->vals, i);
//   for (auto i = 0; i < nmixvals; ++i) Spheral2Silo<T>::readElement(values[i], ucdvar->mixvals, i);

//   // That's it.
//   for (auto i = 0; i != ucdvar->nvals; ++i) {
//     delete[] ucdvar->vals[i];
//     delete[] ucdvar->mixvals[i];
//   }
//   delete[] ucdvar->vals;
//   delete[] ucdvar->mixvals;
//   delete[] ucdvar;
// }

//------------------------------------------------------------------------------
// DBPutUcdvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar1(DBfile& file,
             std::string name,
             std::string meshName,
             std::vector<T>& values,
             std::vector<T>& mixValues,
             int centering,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

  void *values_begin = nullptr, *mix_values_begin = nullptr;
  if (not values.empty()) values_begin = (void*) &(*values.begin());
  if (not mixValues.empty()) mix_values_begin = (void*) &(*mixValues.begin());

  return DBPutUcdvar1(&file,
                      name.c_str(),
                      meshName.c_str(),
                      values_begin,
                      values.size(),
                      mix_values_begin,
                      mixValues.size(),
                      SiloTraits<T>::datatype(),
                      centering,
                      optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBPutQuadvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutQuadvar(DBfile& file,
             std::string name,
             std::string meshName,
             std::vector<T>& values,
             std::vector<T>& mixValues,
             int centering,
             std::vector<int>& vardims,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_ZONECENT);
  auto ndims = vardims.size();
  VERIFY(ndims == 1 or ndims == 2 or ndims == 3);

  // Build the sub-variable names.
  auto nvars = Spheral2Silo<T>::numElements();
  auto nels = values.size();
  auto mixlen = mixValues.size();
  std::vector<char*> varnames;
  for (auto i = 0; i != (int)nvars; ++i) {
    varnames.push_back(const_cast<char*>((name + "_").c_str()));
    sprintf(varnames.back(), "%i", i);
  }

  // Build the sub-variables.
  double** vars = new double*[nvars];
  double** mixvars = new double*[nvars];
  for (auto i = 0; i != (int)nvars; ++i) {
    vars[i] = new double[nels];
    mixvars[i] = new double[mixlen];
  }
  for (auto j = 0; j != (int)nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);
  for (auto j = 0; j != (int)mixlen; ++j) Spheral2Silo<T>::copyElement(mixValues[j], mixvars, j);

  const auto result = DBPutQuadvar(&file,
                                   name.c_str(),
                                   meshName.c_str(),
                                   nvars,
                                   &varnames.front(),
                                   (void*) vars,
                                   &vardims.front(),
                                   ndims,
                                   (void*) mixvars,
                                   mixlen,
                                   SiloTraits<double>::datatype(),
                                   centering,
                                   optlist.mOptlistPtr);

  // That's it.
  for (auto i = 0; i != (int)nvars; ++i) {
    delete[] vars[i];
    delete[] mixvars[i];
  }
  delete[] vars;
  delete[] mixvars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutQuadvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutQuadvar1(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<T>& values,
              std::vector<T>& mixValues,
              int centering,
              std::vector<int>& vardims,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);
  const auto ndims = vardims.size();
  VERIFY(ndims == 1 or ndims == 2 or ndims == 3);

  return DBPutQuadvar1(&file,
                       name.c_str(),
                       meshName.c_str(),
                       (void*) &values.front(),
                       &vardims.front(),
                       ndims,
                       (void*) &mixValues.front(),
                       mixValues.size(),
                       SiloTraits<T>::datatype(),
                       centering,
                       optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBPutZonelist2
//------------------------------------------------------------------------------
inline
int
DBPutZonelist2(DBfile& file,
               std::string name,
               unsigned ndims,
               std::vector<std::vector<int> >& zoneNodes,
               unsigned low_offset,
               unsigned high_offset,
               std::vector<int>& shapetype,
               std::vector<int>& shapesize,
               std::vector<int>& shapecount,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nzones = zoneNodes.size();
  const unsigned nshapes = shapetype.size();
  VERIFY2(shapetype.size() <= nzones, "Bad size: " << shapetype.size() << " !<= " << nzones);
  VERIFY2(shapesize.size() <= nzones, "Bad size: " << shapesize.size() << " !<= " << nzones);
  VERIFY2(shapecount.size() <= nzones, "Bad size: " << shapecount.size() << " !<= " << nzones);
  VERIFY2(shapesize.size() == nshapes, "Bad size: " << shapesize.size() << " != " << nshapes);
  VERIFY2(shapecount.size() == nshapes, "Bad size: " << shapecount.size() << " != " << nshapes);

  // Construct the flat array of zone nodes.
  std::vector<int> nodelist;
  for (unsigned k = 0; k != zoneNodes.size(); ++k) {
    //nodelist.push_back(zoneNodes[k].size());
    std::copy(zoneNodes[k].begin(), zoneNodes[k].end(), std::back_inserter(nodelist));
  }

  // Do the deed.
  const int result = DBPutZonelist2(&file, 
                                    name.c_str(),
                                    nzones,
                                    ndims,
                                    &nodelist.front(),
                                    nodelist.size(),
                                    0,
                                    low_offset,
                                    high_offset,
                                    &shapetype.front(),
                                    &shapesize.front(),
                                    &shapecount.front(),
                                    nshapes,
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutPHZonelist
//------------------------------------------------------------------------------
inline
int
DBPutPHZonelist(DBfile& file,
                std::string name,
                std::vector<std::vector<int> >& faceNodeLists,
                std::vector<std::vector<int> >& zoneFaceLists,
                unsigned low_offset,
                unsigned high_offset,
                DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nfaces = faceNodeLists.size();
  const unsigned nzones = zoneFaceLists.size();

  // Construct the flat arrays of face-node info and zone-face info.
  std::vector<int> nodecnts, nodelist, facecnts, facelist;
  for (unsigned k = 0; k != nfaces; ++k) {
    const std::vector<int>& faceNodes = faceNodeLists[k];
    nodecnts.push_back(faceNodes.size());
    std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(nodelist));
  }
  for (unsigned k = 0; k != nzones; ++k) {
    const std::vector<int>& zoneFaces = zoneFaceLists[k];
    facecnts.push_back(zoneFaces.size());
    std::copy(zoneFaces.begin(), zoneFaces.end(), std::back_inserter(facelist));
  }
  CHECK(nodecnts.size() == nfaces);
  CHECK(facecnts.size() == nzones);

  // Do the deed.
  const int result = DBPutPHZonelist(&file, 
                                     name.c_str(),
                                     nfaces,
                                     &nodecnts.front(),
                                     nodelist.size(),
                                     &nodelist.front(),
                                     NULL,
                                     nzones,
                                     &facecnts.front(),
                                     facelist.size(),
                                     &facelist.front(),
                                     0,
                                     low_offset,
                                     high_offset,
                                     optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointmesh
//------------------------------------------------------------------------------
inline
int
DBPutPointmesh(DBfile& file,
               std::string name,
               std::vector<std::vector<double> >& coords,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned npoints = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == npoints);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[npoints];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Do the deed.
  const int result = DBPutPointmesh(&file, 
                                    name.c_str(),
                                    ndims,
                                    coordPtrs,
                                    npoints,
                                    SiloTraits<double>::datatype(),
                                    optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<T>& values,
              DBoptlist_wrapper& optlist) {

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();

  // Build the sub-variables.
  double** vars = new double*[nvars];
  for (i = 0; (int)i != nvars; ++i) {
    vars[i] = new double[nels];
  }
  for (j = 0; (int)j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);

  const int result = DBPutPointvar(&file,
                                   name.c_str(),
                                   meshName.c_str(),
                                   nvars,
                                   (void*) vars,
                                   nels,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  for (i = 0; (int)i != nvars; ++i) {
    delete[] vars[i];
  }
  delete[] vars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar1(DBfile& file,
               std::string name,
               std::string meshName,
               std::vector<T>& values,
               DBoptlist_wrapper& optlist) {

  return DBPutPointvar1(&file,
                        name.c_str(),
                        meshName.c_str(),
                        (void*) &(*values.begin()),
                        values.size(),
                        SiloTraits<T>::datatype(),
                        optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBAddRegion
//------------------------------------------------------------------------------
inline
int
DBAddRegion(DBmrgtree_wrapper& tree,
            std::string reg_name,
            int info_bits,
            int max_children,
            std::string maps_name,
            std::vector<int>& seg_ids,
            std::vector<int>& seg_lens,
            std::vector<int>& seg_types,
            DBoptlist_wrapper& optlist) {
  int nsegs = seg_ids.size();
  VERIFY((int)seg_lens.size() == nsegs);
  VERIFY((int)seg_types.size() == nsegs);
  return DBAddRegion(tree.mDBmrgtree,
                     reg_name.c_str(),
                     info_bits,
                     max_children,
                     (maps_name.size() > 0 ? maps_name.c_str() : NULL),
                     nsegs,
                     (nsegs > 0 ? &(*seg_ids.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_lens.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_types.begin()) : NULL),
                     optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBSetCwr
//------------------------------------------------------------------------------
inline
int
DBSetCwr(DBmrgtree_wrapper& tree,
         std::string path) {
  return DBSetCwr(tree.mDBmrgtree,
                  path.c_str());
}

//------------------------------------------------------------------------------
// DBGetCwr
//------------------------------------------------------------------------------
inline
const char*
DBGetCwr(DBmrgtree_wrapper& tree) {
  return DBGetCwr(tree.mDBmrgtree);
}

//------------------------------------------------------------------------------
// DBPutMrgtree
//------------------------------------------------------------------------------
inline
int
DBPutMrgtree(DBfile& file,
             std::string name,
             std::string mesh_name,
             DBmrgtree_wrapper& tree,
             DBoptlist_wrapper& optlist) {
  return DBPutMrgtree(&file,
                      name.c_str(),
                      mesh_name.c_str(),
                      tree.mDBmrgtree,
                      optlist.mOptlistPtr);
}

}

//------------------------------------------------------------------------------
// Typedef for vectors of optlists.
//------------------------------------------------------------------------------
typedef std::vector<silo::DBoptlist_wrapper*> vector_of_DBoptlist;

#endif
