//
//  EulerianStencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef EulerianStencil_h
#define EulerianStencil_h

#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/LinearAlgebra/SquareCube.h"
#include "JMM_CPPLibs/Macros/TemplateLog2.h"
#include "CommonStencil.h"
// ------------ Eulerian scheme -------------

template<typename TS, size_t n> struct QuadraticMax;

// --- Stencil --
/** A PDE discretization is specified via a set of differences, of several types,
 which number is fixed in advance by the following parameters.*/

template<typename TDiff, int nSym, int nFor=0, int nM=1> struct
EulerianStencil {
    typedef TDiff DifferenceType;
    static const int nForward = nFor, nSymmetric = nSym, nMax = nM;
    Redeclare3Types(DifferenceType, ScalarType,OffsetType,MultiplierType)

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
	struct ActiveNeighFlagType : std::bitset<nNeigh+nMaxBits> {
		using Superclass = std::bitset<nNeigh+nMaxBits>;
		explicit operator ScalarType() const {return this->to_ullong();}
		explicit ActiveNeighFlagType(ScalarType a):Superclass((unsigned long long)a){}
		ActiveNeighFlagType(){}
	};
    static int GetIMax(ActiveNeighFlagType b);
    static void SetIMax(ActiveNeighFlagType & b, int iMax);
    
    typedef QuadraticMax<ScalarType,nMax> QuadType;
    ScalarType HopfLaxUpdate(OffsetType, ScalarType, const MultiplierType &,
                             QuadType &, ActiveNeighFlagType &) const;
	
	using CommonStencilType = CommonStencil<OffsetType,ScalarType,nActiveNeigh>;
	Redeclare3Types(CommonStencilType,DiscreteFlowElement,DiscreteFlowType,RecomputeType);
    template<typename F> RecomputeType
	HopfLaxRecompute(const F &,const MultiplierType &,
					 ActiveNeighFlagType, DiscreteFlowType &) const;
};

// ---- Enhanced offsets, referred to as differences  ---
/** A Difference is a basic component of a PDE scheme. It is the data of an offset and weight.
 The weight which is either specified directly or as a baseweight and a multiplier index, within [0,VMultSize[.*/

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
    ScalarType Weight(const MultiplierType & mult) const {
		return baseWeight*square(mult[multIndex]);}
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
    ScalarType Weight(const MultiplierType & mult) const {
		return baseWeight*square(mult);}
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
    static constexpr ScalarType Infinity() {
		return std::numeric_limits<ScalarType>::infinity();}
    
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


#include "Implementation/EulerianStencil.hxx"

#endif /* EulerianStencil_h */
