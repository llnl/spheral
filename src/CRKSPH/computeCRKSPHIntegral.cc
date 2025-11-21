//------------------------------------------------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#include <stdio.h>
#include "computeCRKSPHEvaluation.hh"
#include "computeCRKSPHIntegral.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
std::pair<typename Dimension::Vector,typename Dimension::Vector>
computeCRKSPHIntegral(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& weight,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       size_t nodeListi, const int i, size_t nodeListj, const int j, int mydim, const int order,
                       typename Dimension::Vector rmin, typename Dimension::Vector rmax){


  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  bool is_zero=false;//flag used to exit routine if integral is zero, to increase efficiency 

  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Pre-conditions.
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  const Scalar wi = weight(nodeListi, i);
  const Scalar wj = weight(nodeListj, j);

  //Vector dumvec=Vector::zero();//Dummy Variables
  //Scalar dumscal=0.0;
  Vector temprmin=rmin;
  Vector temprmax=rmax;

  int dim = Dimension::nDim;
  //Vector integ[order][order];
  //Vector integ2[order][order];
  vector< vector<Vector> > integ(order, vector<Vector>(order));
  vector< vector<Vector> > integ2(order, vector<Vector>(order));
  Scalar Wai=0.0; 
  Scalar Waj=0.0; 
  Vector gradWai=Vector::zero(); 
  Vector gradWaj=Vector::zero(); 
  if(dim == mydim){//Lowest level perform 1d normal Romberg

    computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListi, i, rmin,
                          true, Wai, gradWai);
    computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListj, j, rmin,
                          true, Waj, gradWaj);
 
    integ[0][0]=Wai*gradWaj*wi*wj;
    integ2[0][0]=Waj*gradWai*wi*wj;
    computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListi, i, rmax,
                        true, Wai, gradWai);
    computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListj, j, rmax,
                        true, Waj, gradWaj);
    integ[0][0]+=Wai*gradWaj*wi*wj;
    integ2[0][0]+=Waj*gradWai*wi*wj;
  }else{
    temprmin=rmin;
    temprmax=rmax;
    temprmax(mydim-1) = rmin(mydim-1);
    integ[0][0]=computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).first;
    integ2[0][0]=computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).second;
    temprmax=rmax;
    temprmin(mydim-1) = rmax(mydim-1);
    integ[0][0]+=computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).first;
    integ2[0][0]+=computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).second;
  }
  integ[0][0]*=(rmax(mydim-1)-rmin(mydim-1))*0.5;
  integ2[0][0]*=(rmax(mydim-1)-rmin(mydim-1))*0.5;
  
  for(int nn=1; nn < order; nn++){
    for(int mm=0; mm <= nn; mm++){
       if(mm==0){
          integ[nn][mm]=integ[nn-1][mm]*0.5;
          integ2[nn][mm]=integ[nn-1][mm]*0.5;
          Scalar hn = (rmax(mydim-1)-rmin(mydim-1))/pow(2.0,nn);
          for(int k=1; k<= pow(2.0,nn-1); k++){
              Vector reval=rmin;
              reval(mydim-1)+=(2*k-1)*hn;
              if(dim == mydim){
                 computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListi, i, reval,
                           true, Wai, gradWai);
                 computeCRKSPHEvaluation(connectivityMap,W,weight,position,H,nodeListj, j, reval,
                           true, Waj, gradWaj);
                 integ[nn][mm]+=hn*Wai*gradWaj*wi*wj;
                 integ2[nn][mm]+=hn*Waj*gradWai*wi*wj;
              }else{
                temprmin=rmin;
                temprmax=rmax;
                temprmax(mydim-1)=reval(mydim-1);
                temprmin(mydim-1)=reval(mydim-1);
                integ[nn][mm]+=hn*computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).first;
                integ2[nn][mm]+=hn*computeCRKSPHIntegral(connectivityMap,W,weight,position,H,nodeListi, i, nodeListj, j, mydim+1, order, temprmin, temprmax).second;
              }
             
          }
       }else{
         integ[nn][mm]=(1.0/(pow(4.0,mm)-1.0))*(pow(4.0,mm)*integ[nn][mm-1]-integ[nn-1][mm-1]);
         integ2[nn][mm]=(1.0/(pow(4.0,mm)-1.0))*(pow(4.0,mm)*integ[nn][mm-1]-integ[nn-1][mm-1]);
       }
 
    }
  }
  return std::pair<Vector,Vector>(integ[order-1][order-1],integ2[order-1][order-1]);

}

}

