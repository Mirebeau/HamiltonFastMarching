//
//  StencilType.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef StencilType_h
#define StencilType_h

#include "Output/ExportMacros.h"

// ------------ Eulerian scheme -------------

// Silly replacement for a constexpr log. TODO : hide somewhere.

template<size_t i> struct FloorLog2 {static const size_t value = 1+FloorLog2<i/2>::value;};
template<> struct FloorLog2<1> {static const size_t value = 0;};
template<size_t i> class CeilLog2 {
    static const size_t floor = FloorLog2<i>::value;
public:
    static const size_t value = (i== (1<<floor) ? floor : (floor+1));
};


// Stencil
template<typename TDiff, int nSym, int nFor=0, int nM=1> struct
EulerianStencil {
    typedef TDiff DifferenceType;
    static const int nForward = nFor, nSymmetric = nSym, nMax = nM;

    std::array<std::array<DifferenceType, nSymmetric>, nMax> symmetric;
    std::array<std::array<DifferenceType, nForward>, nMax> forward;
    PrintSelfMacro(EulerianStencil);

    /// Largest number of simultaneously active neighbors
    static const int nActiveNeigh = nForward+nSymmetric;
    /// Number of neighbors for a single value of iMax
    static const int nSingleNeigh = nForward+2*nSymmetric;
    /// Total number of neighbors
    static const int nNeigh = nMax*nSingleNeigh;
    
    /// Flat to store simultaneously active neighbors
    static const int nMaxBits = CeilLog2<nMax>::value;
    typedef std::bitset<nNeigh+nMaxBits> ActiveNeighFlagType;
    static int GetIMax(ActiveNeighFlagType b) {
        int iMax=0;
        for(int i=0; i<nMaxBits; ++i) if(b[nNeigh+i]) iMax += (1<<i);
        assert(0<=iMax && iMax<nMax);
        return iMax;
    }
    static void SetIMax(ActiveNeighFlagType & b, int iMax) {
        assert(0<=iMax && iMax<nMax);
        for(int i=0; i<nMaxBits; ++i) b[nNeigh+i] = iMax & (1<<i);
    }
};

// Enhanced offsets, referred to as differences (define the finite difference scheme)
template<typename TOff, typename TScal, int VMult>
struct EulerianDifference {
    typedef TOff OffsetType;
    typedef typename OffsetType::ComponentType ShortType;

    typedef TScal ScalarType;
    static const int multSize=VMult;
    typedef LinearAlgebra::Point<ScalarType, multSize> MultiplierType;
    
    OffsetType offset;
    ScalarType baseWeight;
    ShortType multIndex;
    ScalarType Weight(const MultiplierType & mult) const {return baseWeight*square(mult[multIndex]);}
};

template<typename TOff, typename TScal>
struct EulerianDifference<TOff,TScal,1> {
    typedef TOff OffsetType;
    typedef typename OffsetType::ComponentType ShortType;
    
    typedef TScal ScalarType;
    static const int multSize=1;
    typedef ScalarType MultiplierType;

    OffsetType offset;
    ScalarType baseWeight;
    ScalarType Weight(const MultiplierType & mult) const {return baseWeight*square(mult);}
};

template<typename TOff, typename TScal>
struct EulerianDifference<TOff,TScal,0> {
    typedef TOff OffsetType;
    typedef typename OffsetType::ComponentType ShortType;
    
    typedef TScal ScalarType;
    static const int multSize=0;
    struct MultiplierType {};
    
    OffsetType offset;
    ScalarType baseWeight;
    ScalarType Weight(MultiplierType) const {return baseWeight;}
};

// Printing


template<typename TDiff, int nFor, int nSym, int nM> void
EulerianStencil<TDiff,nFor,nSym,nM>::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayRecursiveArrow(forward, 1)
    ExportArrayRecursiveArrow(symmetric, 1)
    << "}";
}

template<typename TOff, typename TScal, int VMult> std::ostream &
operator << (std::ostream & os, const EulerianDifference<TOff, TScal, VMult> & diff){
    return os << "{" << diff.baseWeight << "," << diff.offset  << "}";
}

// ----------- Semi-Lagrangian scheme ------------

enum class LagrangianStencilPeriodicity {None, Simple, Double};

template<typename TOff, typename TIndDiff, LagrangianStencilPeriodicity VPer>
struct LagrangianStencil {
    typedef TOff OffsetType;
    typedef TIndDiff IndexDiff;
    typedef typename IndexDiff::ComponentType DiscreteType;
    static const DiscreteType Dimension = IndexDiff::Dimension;
    static const LagrangianStencilPeriodicity Periodicity;
    
    DiscreteType NSectors(){
        switch(Periodicity){
            case LagrangianStencilPeriodicity::None:    return nOffsets-1;
            case LagrangianStencilPeriodicity::Simple:  return nOffsets;
            case LagrangianStencilPeriodicity::Double:  return 2*nOffsets;
        }
    }
    IndexDiff Sector(DiscreteType n, DiscreteType k){
        assert(0<=n && n<=NSectors())
        assert(0<=k && k<Dimension);
        const DiscreteType m = n+k;
        switch(Periodicity){
            case LagrangianStencilPeriodicity::None: return IndexDiff::CastCoordinates(_begin[m]);
            case LagrangianStencilPeriodicity::Simple: return IndexDiff::CastCoordinates(_begin[m%nOffsets]);
            case LagrangianStencilPeriodicity::Double: return IndexDiff::CastCoordinates(m<nOffsets ? _begin[m] : -_begin[m-nOffsets]);
        }
    }
    
protected:
    OffsetType * _begin;
    DiscreteType nOffsets;
};

#endif /* StencilType_h */
