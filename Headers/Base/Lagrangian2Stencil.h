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
#include "JMM_CPPLibs/LinearAlgebra/RanderNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/AsymmetricQuadraticNorm.h"

#include "Specializations/CommonTraits.h"
#include "JMM_CPPLibs/LinearAlgebra/HopfLaxMinimize.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"

// ----------- Semi-Lagrangian scheme ------------

template<typename TOff, typename TIndDiff>
struct Lagrangian2Stencil {
    typedef TOff OffsetType;
    typedef TIndDiff IndexDiff;
    typedef typename IndexDiff::ComponentType DiscreteType;
    typedef typename OffsetType::ComponentType ShortType;
    static const DiscreteType Dimension = IndexDiff::Dimension;
    
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
    typedef double ScalarType;
    
    struct DiscreteFlowElement {OffsetType offset; ScalarType weight;};
    typedef CappedVector<DiscreteFlowElement, nActiveNeigh> DiscreteFlowType;
    struct RecomputeType {ScalarType value,width;};

//protected:
    OffsetType * pOffsets;
    DiscreteType nOffsets;
};


template<typename TO, typename TID> void
Lagrangian2Stencil<TO,TID>::PrintSelf(std::ostream & os) const {
    os << "{";
    for(int i=0; i<nOffsets; ++i) os << pOffsets[i] << ",";
    os << "}";
}

// -----

struct TraitsRanderLag2 : TraitsBase<2> {
    typedef Lagrangian2Stencil<OffsetType,IndexDiff> StencilType;
    typedef PeriodicGrid<TraitsRanderLag2> DomainType;
    struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
    
    typedef LinearAlgebra::RanderNorm<ScalarType,2> NormType;
    typedef LinearAlgebra::RanderNorm<ScalarType,1> NormType1;
    typedef NormType DistanceGuess;
};

struct TraitsAsymmetricQuadraticLag2 : TraitsBase<2> {
    typedef Lagrangian2Stencil<OffsetType,IndexDiff> StencilType;
    typedef PeriodicGrid<TraitsAsymmetricQuadraticLag2> DomainType;
    struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
    
    typedef LinearAlgebra::AsymmetricQuadraticNorm<ScalarType,2> NormType;
    typedef LinearAlgebra::AsymmetricQuadraticNorm<ScalarType,1> NormType1;
    typedef NormType DistanceGuess;
};

template<typename TPred, typename TVec>
void SternBrocotRefine(const TPred & stop, std::forward_list<TVec> & l){
    typedef TVec VectorType;
    assert(!l.empty());
    int size = 8; // Typical upper bound. TODO : actual size ?
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


/// Stencil for metrics built of a quadratic and a linear part
template<typename T>
struct StencilQuadLinLag2 : HamiltonFastMarching<T>::StencilDataType {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::StencilDataType Superclass;
	Redeclare6Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,RecomputeType,Traits)
    Redeclare7Types(HFM,IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,DomainType,IndexDiff)

Redeclare6Types(Traits,NormType,NormType1,IndexType,StencilType,OffsetType,DistanceGuess)
	
	Redeclare1Type(Superclass,OffsetVal3)
    Redeclare1Constant(HFM,Dimension)
    
    // Specific to this model
    virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVal3 &) override;
    virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;

    typedef typename NormType::SymmetricMatrixType SymmetricMatrixType;
    typedef LinearAlgebra::VectorPair<SymmetricMatrixType,VectorType> MetricElementType; // Allows averaging
    
    // Generic
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pMetric;
    ParamDefault param;

    virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *) override;
    virtual DistanceGuess GetGuess(IndexCRef index) const override;
private:
    std::forward_list<OffsetType> l;
    NormType GetNorm(IndexCRef index) const {
        const MetricElementType data = (*pMetric)(index);
        return NormType{data.first,data.second};}
    
    // Refinement near walls
    ScalarType wallBoundaryAngularResolution = 0.2;
    typename Traits::template Array<bool,Dimension> walls;
/*    typedef typename HFM::template DataSource<bool> BoolField;
    std::unique_ptr<BoolField> pWalls = nullptr;*/
    const DomainType * pDom = nullptr;
    bool OnWallBoundary(IndexCRef) const;
};


typedef StencilQuadLinLag2<TraitsRanderLag2> StencilRanderLag2;
typedef StencilQuadLinLag2<TraitsAsymmetricQuadraticLag2> StencilAsymmetricQuadraticLag2;

#include "Lagrangian2Stencil.hxx"

#endif /* Lagrangian2Stencil_h */
