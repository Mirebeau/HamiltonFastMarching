//
//  Lagrangian2Stencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef Lagrangian2Stencil_h
#define Lagrangian2Stencil_h

#include <forward_list>
#include "JMM_CPPLibs/DataStructures/CappedVector.h"

// ----------- Semi-Lagrangian scheme ------------

template<typename TOff, typename TScalar, typename TDiscrete>
struct Lagrangian2Stencil {
    typedef TOff OffsetType;
	typedef TDiscrete DiscreteType;
	typedef TScalar ScalarType;
    typedef typename OffsetType::ComponentType ShortType;
    static const DiscreteType Dimension = OffsetType::Dimension;
	static_assert(Dimension==2,"Two dimensional stencil class");
    
	DiscreteType NSectors() const {return nOffsets;}
    const OffsetType & Sector(DiscreteType n, DiscreteType k) const {
        assert(0<=k && k<Dimension);
        assert(0<=n && n<NSectors()); // slightly abusive use case
        return pOffsets[(n+k)%nOffsets];
    }
    
    PrintSelfMacro(Lagrangian2Stencil);
    struct ActiveNeighFlagType {
        ShortType sectorIndex;
        ActiveNeighFlagType():sectorIndex(-1){};
        ActiveNeighFlagType(long a):sectorIndex(ShortType(a)){};
        bool none() const {return sectorIndex==-1;}
        unsigned long to_ulong() const {return (unsigned long)(sectorIndex);}
    };
    
    // TODO : Remove when possible ? (Put in StencilData specialization ?)
    static const int nActiveNeigh = Dimension;
    struct DiscreteFlowElement {OffsetType offset; ScalarType weight;};
    typedef CappedVector<DiscreteFlowElement, nActiveNeigh> DiscreteFlowType;
    struct RecomputeType {ScalarType value,width;};

    OffsetType * pOffsets;
    DiscreteType nOffsets;
};


template<typename TO, typename TS, typename TD> void
Lagrangian2Stencil<TO,TS,TD>::PrintSelf(std::ostream & os) const {
    os << "{";
    for(int i=0; i<nOffsets; ++i) os << pOffsets[i] << ",";
    os << "}";
}

// -----

// TODO : use two vectors instead of a forward list ?
template<typename TPred, typename TVec>
void SternBrocotRefine(const TPred & stop, std::forward_list<TVec> & l){
    typedef TVec VectorType;
    assert(!l.empty());
    int size = 8; // There is no l.size(), but this is a typical upper bound.
    auto it0 = l.begin(), it1=it0; ++it1;
    for(; it1!=l.end(); ){
        const VectorType & u=*it0, & v=*it1;
        if(stop(u,v)){
            ++it0;
            ++it1;
        } else {
            it1 = l.insert_after(it0,u+v);
            ++size;
            if(size>=127){ExceptionMacro("Stern-Brocot refine error : excessive stencil size. "
                                         "Metric is either non-definite or too anisotropic.");}
        }            
    }
}

#endif /* Lagrangian2Stencil_h */
