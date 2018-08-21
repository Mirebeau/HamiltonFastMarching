//
//  EulerianStencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef EulerianStencil_h
#define EulerianStencil_h

#include "Output/ExportMacros.h"
#include "LinearAlgebra/SquareCube.h"

// ------------ Eulerian scheme -------------

template<typename TS, size_t n> struct QuadraticMax;

// Silly replacement for a constexpr log. TODO : hide somewhere.

template<size_t i> struct FloorLog2 {static const size_t value = 1+FloorLog2<i/2>::value;};
template<> struct FloorLog2<1> {static const size_t value = 0;};
template<size_t i> class CeilLog2 {
    static const size_t floor = FloorLog2<i>::value;
public:
    static const size_t value = (i== (1<<floor) ? floor : (floor+1));
};


// --- Stencil --

#define FromDifferenceType(x) DifferenceType:: x

template<typename TDiff, int nSym, int nFor=0, int nM=1> struct
EulerianStencil {
    typedef TDiff DifferenceType;
    static const int nForward = nFor, nSymmetric = nSym, nMax = nM;

    // All data
    std::array<std::array<DifferenceType, nSymmetric>, nMax> symmetric;
    std::array<std::array<DifferenceType, nForward>, nMax> forward;
    PrintSelfMacro(EulerianStencil);

    
    /// Largest number of simultaneously active neighbors
    static const int nActiveNeigh = nForward+nSymmetric;
    /// Number of neighbors for a single value of iMax
    static const int nSingleNeigh = nForward+2*nSymmetric;
    /// Total number of neighbors
    static const int nNeigh = nMax*nSingleNeigh;
    static const int nMaxBits = CeilLog2<nMax>::value;
    
    /// Encodes which neighbors are active at a given point
    typedef std::bitset<nNeigh+nMaxBits> ActiveNeighFlagType;
    static int GetIMax(ActiveNeighFlagType b);
    static void SetIMax(ActiveNeighFlagType & b, int iMax);
    
    Redeclare3Types(FromDifferenceType, ScalarType,OffsetType,MultiplierType)
    typedef QuadraticMax<ScalarType,nMax> QuadType;
    ScalarType HopfLaxUpdate(OffsetType, ScalarType, const MultiplierType &,
                             QuadType &, ActiveNeighFlagType &) const;
    
    struct DiscreteFlowElement {OffsetType offset; ScalarType weight;};
    typedef CappedVector<DiscreteFlowElement, nActiveNeigh> DiscreteFlowType;
    struct RecomputeType {ScalarType value,width;};
    template<typename F> RecomputeType HopfLaxRecompute(const F &, const MultiplierType &, ActiveNeighFlagType, DiscreteFlowType &) const;
};

// ---- Enhanced offsets, referred to as differences (define the finite difference scheme) ---
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


// ---- Data structure related to the max of a series of quadratic form. -----
// ----- Used to solve the numerical scheme -----

template<typename TS, size_t n>
struct QuadraticMax {
    typedef TS ScalarType;
    static constexpr ScalarType Infinity() {return std::numeric_limits<ScalarType>::infinity();}
    
    ScalarType minVal=-Infinity();
    
    std::pair<ScalarType, int> Solve() const;
    void Add(ScalarType value, ScalarType weight, int);
    void Add(ScalarType value, ScalarType weight);
    PrintSelfMacro(QuadraticMax);
protected:
    // Note : result is not really needed, except for some assertions.
    struct DataType {ScalarType a=0, b=0, c=-1; PrintSelfMacro(DataType);};
    std::array<DataType,n> data;
};


#include "EulerianStencil.hxx"

#endif /* EulerianStencil_h */
