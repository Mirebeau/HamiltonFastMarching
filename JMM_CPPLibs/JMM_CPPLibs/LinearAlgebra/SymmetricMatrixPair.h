// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef LiftedFastMarching_SymmetricMatrixPair_h
#define LiftedFastMarching_SymmetricMatrixPair_h

#include <type_traits>
#include "VectorType.h"

namespace LinearAlgebra {
    
    template<
    typename TSymmetricMatrix0,
    typename TSymmetricMatrix1
    >
    struct SymmetricMatrixPair :
    vector_space<
    SymmetricMatrixPair<TSymmetricMatrix0, TSymmetricMatrix1>,
    typename TSymmetricMatrix0::ComponentType
    >
{
    typedef typename TSymmetricMatrix0::ComponentType ComponentType;
    
    // ** Miscellaneous typedefs **
    
    typedef TSymmetricMatrix0 SymmetricMatrix0;
    typedef TSymmetricMatrix1 SymmetricMatrix1;
    
    SymmetricMatrix0 m0;
    SymmetricMatrix1 m1;
    SymmetricMatrixPair(const SymmetricMatrix0 & m0_, const SymmetricMatrix1 & m1_):m0(m0_),m1(m1_){}
    SymmetricMatrixPair(){};

    
    typedef typename SymmetricMatrix0::VectorType VectorType0;
    typedef typename SymmetricMatrix1::VectorType VectorType1;
    
    static const int Dimension0 = VectorType0::Dimension;
    static const int Dimension1 = VectorType1::Dimension;
    static const int Dimension = Dimension0+Dimension1;
    typedef Vector<ComponentType, Dimension> VectorType;
    
    // ** Splitting vectors and points into components **
    // Dummy variables are introduced in the dimensionless case
    
    static       VectorType0 & Part0(      VectorType & u) {static VectorType0 DummyVector;
        return Dimension0==0 ? DummyVector : *reinterpret_cast<      VectorType0*>(&u[0]);}
    static const VectorType0 & Part0(const VectorType & u) {static VectorType0 DummyVector;
        return Dimension0==0 ? DummyVector : *reinterpret_cast<const VectorType0*>(&u[0]);}
    static       VectorType1 & Part1(      VectorType & u) {static VectorType1 DummyVector;
        return Dimension1==0 ? DummyVector : *reinterpret_cast<      VectorType1*>(&u[Dimension0]);}
    static const VectorType1 & Part1(const VectorType & u) {static VectorType1 DummyVector;
        return Dimension1==0 ? DummyVector : *reinterpret_cast<const VectorType1*>(&u[Dimension0]);}
    
    typedef typename SymmetricMatrix0::PointType PointType0;
    typedef typename SymmetricMatrix1::PointType PointType1;
    typedef Point<ComponentType, Dimension> PointType;

    static       PointType0 & Part0(      PointType & u) {static PointType0 DummyPoint;
        return Dimension0==0 ? DummyPoint : *reinterpret_cast<      PointType0*>(&u[0]);}
    static const PointType0 & Part0(const PointType & u) {static PointType0 DummyPoint;
        return Dimension0==0 ? DummyPoint : *reinterpret_cast<const PointType0*>(&u[0]);}
    static       PointType1 & Part1(      PointType & u) {static PointType1 DummyPoint;
        return Dimension1==0 ? DummyPoint : *reinterpret_cast<      PointType1*>(&u[Dimension0]);}
    static const PointType1 & Part1(const PointType & u) {static PointType1 DummyPoint;
        return Dimension1==0 ? DummyPoint : *reinterpret_cast<const PointType1*>(&u[Dimension0]);}

    // ** Block diagonal matrix algebra **
    
    ComponentType ScalarProduct(const VectorType & u, const VectorType & v) const {
        return
        m0.ScalarProduct(Part0(u),Part0(v))+
        m1.ScalarProduct(Part1(u),Part1(v));
    }
    ComponentType SquaredNorm(const VectorType & u) const {return ScalarProduct(u, u);}
    ComponentType Norm(const VectorType & u) const {return sqrt(SquaredNorm(u));}
    
    SymmetricMatrixPair Inverse() const {
        SymmetricMatrixPair result;
        result.m0 = m0.Inverse();
        result.m1 = m1.Inverse();
        return result;
    }
    
    VectorType operator*(const VectorType & u) const {
        VectorType v;
        Part0(v) = m0*Part0(u);
        Part1(v) = m1*Part1(u);
        return v;
    }
	
protected:
    static_assert(std::is_same<ComponentType, typename TSymmetricMatrix1::ComponentType>::value && std::is_same<ComponentType, typename VectorType0::ComponentType>::value && std::is_same<ComponentType, typename VectorType1::ComponentType>::value,
                  "ComponentTypes must match");
};
    
    template<typename TS0, typename TS1>
    std::ostream & operator << (std::ostream & f, const SymmetricMatrixPair<TS0, TS1> & m)
    {
        f << "{" << m.m0 << ", " << m.m1 << "}";
        return f;
    }
}


#endif
