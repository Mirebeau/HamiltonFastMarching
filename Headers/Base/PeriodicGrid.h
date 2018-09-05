// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef PeriodicGrid_h
#define PeriodicGrid_h

#include "JMM_CPPLibs/DataStructures/RedeclareTypesMacro.h"
#include "JMM_CPPLibs/Output/ExceptionMacro.h"
#include "JMM_CPPLibs/LinearAlgebra/SquareCube.h"

// Spherical boundary conditions correspond to the parametrization (here in dim 4)
//( cos(t1), sin(t1) cos(t2), sin(t1) sin(t2) cos(t3), sin(t1) sin(t2) sin(t3))
// which is invariant under the transformations ti->2pi-ti, t(i+1)->t(i+1)+pi.
// with 0<=ti<=Pi except for the last 0<= tn <=2Pi, and potentially the first 0<=t1<=Pi/2 for the projective space.
// Obviously the adequate metric must be used.
// This coordinate system is also suitable for the projective space.

enum class Boundary {
    Closed,     // Bounded box
    Periodic,   // Standard periodicity
    Sphere2_0, Sphere2_1, Sphere2_Hopf, // 2D Sphere (with circle fibration)
    Sphere3_0, Sphere3_1, Sphere3_2, Proj3_1 // 3D Sphere or 3D Projective space
//    Sphere, SphereLast
};

struct Boundary_AllClosed {
    Boundary operator[](size_t) const {return Boundary::Closed;}
};

/* 2D sphere boundary conditions
 Meant for parametrization
 (t0,t1)->(cos(t0), sin(t0)cos(t1),sin(t0)sin(t1))
 with 0<=t0<=pi, 0<=t1<=2pi, n1=2*n0,
 and metric dt0^2+sin(t0)^2 dt1^2.

 Implements the periodicity relations (t0,t1)~(t0,t1+2pi) and (t0,t1)~(2pi-t0,t1+pi).
 
 Can also be used for 2D projective space, with 0<=t0<=pi/2, n1=4n0.
 */

 /* 3D Sphere boundary conditions
  Meant for parametrization
  (t0,t1,t2)->(cos(t0)cos(t1), cos(t0)sin(t1), sin(t0)cos(t2),sin(t0)sin(t2))
  with 0<=t0<=pi/2, 0<=t1<=2pi, 0<=t2<=2pi, n1=n2=4n0.
  and metric dt0^2+cos(t0)^2 dt1^2+sin(t0)^2dt2^2.
  
  Periodicity: (t0,t1,t2)~(t0+pi,t1+pi,t2+pi), (pi-t0,t1+pi,t2), (t0,t1+2pi,t2), (t0,t1,t2+2pi).
  
  Can also be used for 3D projective space, with 0<=t1<=pi, n1=2n0, (t0,t1,t2)~(t0,t1+pi,t2+pi).
  */

template<typename TTraits> struct PeriodicGrid {
    typedef TTraits Traits;
    Redeclare4Types(FromTraits,DiscreteType,ScalarType,IndexType,PointType);
    Redeclare1Constant(FromTraits,Dimension);
    typedef const IndexType & IndexCRef;
    
    DiscreteType LinearFromIndex(IndexCRef) const;
    IndexType    IndexFromLinear(DiscreteType) const;
    PointType    PointFromIndex(IndexCRef) const;
    IndexType    IndexFromPoint(const PointType &) const;
    
    // Sets back to parametrization domain and returns which axes were flipped + failbit.
    typedef std::bitset<Dimension+1> ReverseFlag;
    static constexpr bool MayReverse(DiscreteType i);
    ReverseFlag Periodize(IndexType &) const;
    ReverseFlag Periodize(PointType &) const;
    
    PeriodicGrid(IndexCRef);
protected:
    typename Traits::template Array<ScalarType,Dimension> arr;
};

// A second class is used to reparametrize the domain
template<typename TPoint, typename TVec> struct ParamInterface_ {
    typedef TPoint PointType;
    typedef TVec VectorType;
    virtual PointType ADim(const PointType & p) const=0;
    virtual PointType ReDim(const PointType & p) const=0;
    virtual VectorType ReDim(const VectorType & p) const=0;
    virtual ~ParamInterface_(){};
};



#include "PeriodicGrid.hxx"

#endif /* PeriodicGrid_h */
