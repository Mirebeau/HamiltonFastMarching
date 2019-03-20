// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef PeriodicGrid_h
#define PeriodicGrid_h

#include "JMM_CPPLibs/Macros/RedeclareTypes.h"
#include "JMM_CPPLibs/Macros/Exception.h"
#include "JMM_CPPLibs/LinearAlgebra/SquareCube.h"

// TODO : Put this file in specializations

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
    constexpr Boundary operator[](int) const {return Boundary::Closed;}
};

template<int VDim, Boundary VLast>
struct Boundary_ClosedButLast {
	constexpr static const int Dimension = VDim;
	constexpr static const Boundary lastBoundary = VLast;
	constexpr Boundary operator[](int i) const {
		assert(0<=i && i<Dimension);
		return i<Dimension-1 ? Boundary::Closed : lastBoundary;
	}
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

template<typename TTraits> struct PeriodicGrid : TTraits::BaseDomain {
    typedef TTraits Traits;
    typedef typename Traits::BaseDomain Superclass;
	Redeclare6Types(Superclass,IndexType,PointType,IndexCRef,DiscreteType,ScalarType,NeighborsType);
    Redeclare1Constant(Superclass,Dimension);    
    struct Transform;
    
    // Target point, base point (here ignored)
    static const bool periodizeUsesBase = false;
    Transform Periodize(IndexType & target, IndexCRef base) const;
    Transform Periodize(PointType & target, const PointType & base) const;
    Transform PeriodizeNoBase(IndexType & target) const;
    Transform PeriodizeNoBase(PointType & target) const;
	
	NeighborsType Neighbors(const PointType &, bool=false) const; // Neighbors on the cartesian grid
	
    PeriodicGrid(IndexCRef);
};


template<typename TTraits> struct
PeriodicGrid<TTraits>::Transform {
    bool IsValid() const {return !reverseFlag[Dimension];}
	bool IsTrivial() const {return reverseFlag.none();}
    template<typename TVec> void PullVector(TVec &) const;
    void Invalidate(){reverseFlag[Dimension]=true;}
protected:
    friend PeriodicGrid<TTraits>;
    std::bitset<Dimension+1> reverseFlag;
    static constexpr bool MayReverse(DiscreteType i);
};

#include "Implementation/PeriodicGrid.hxx"

#endif /* PeriodicGrid_h */
