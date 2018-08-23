//
//  Lagrangian2Stencil.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef Lagrangian2Stencil_h
#define Lagrangian2Stencil_h

#include <forward_list>
#include "DataStructures/CappedVector.h"
#include "LinearAlgebra/RanderNorm.h"
#include "LinearAlgebra/AsymmetricQuadraticNorm.h"

#include "Specializations/CommonTraits.h"
#include "LinearAlgebra/HopfLaxMinimize.h"
#include "LinearAlgebra/VectorPairType.h"

// ----------- Semi-Lagrangian scheme ------------

//enum class Lagrangian2StencilPeriodicity {None, Simple, Double};

template<typename TOff, typename TIndDiff, Lagrangian2StencilPeriodicity VPer>
struct Lagrangian2Stencil {
    typedef TOff OffsetType;
    typedef TIndDiff IndexDiff;
    typedef typename IndexDiff::ComponentType DiscreteType;
    typedef typename OffsetType::ComponentType ShortType;
    static const DiscreteType Dimension = IndexDiff::Dimension;
    static const Lagrangian2StencilPeriodicity Periodicity = VPer;
    
    DiscreteType NSectors() const {
        switch(Periodicity){
            case Lagrangian2StencilPeriodicity::None:    return nOffsets-1;
            case Lagrangian2StencilPeriodicity::Simple:  return nOffsets;
            case Lagrangian2StencilPeriodicity::Double:  return 2*nOffsets;
        }
    }
    OffsetType Sector(DiscreteType n, DiscreteType k) const {
        assert(0<=k && k<Dimension);
        assert(0<=n && n<NSectors()
               || (Periodicity==Lagrangian2StencilPeriodicity::None && n==NSectors() && k==0)); // slightly abusive use case
        const DiscreteType m = n+k;
        switch(Periodicity){
            case Lagrangian2StencilPeriodicity::None:   return pOffsets[m];
            case Lagrangian2StencilPeriodicity::Simple: return pOffsets[m%nOffsets];
            case Lagrangian2StencilPeriodicity::Double: return m<nOffsets ? pOffsets[m] : -pOffsets[m-nOffsets];
        }
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


//typedef Lagrangian2StencilPeriodicity LSP;
template<typename TO, typename TID, Lagrangian2StencilPeriodicity VPer> void
Lagrangian2Stencil<TO,TID,VPer>::PrintSelf(std::ostream & os) const {
    os << "{";
    for(int i=0; i<nOffsets; ++i) os << pOffsets[i] << ",";
    os << "}";
}

// -----

struct TraitsRanderLag2 : TraitsBase<2> {
    typedef Lagrangian2Stencil<OffsetType,IndexDiff,Lagrangian2StencilPeriodicity::Simple> StencilType;
    typedef PeriodicGrid<TraitsRanderLag2> DomainType;
    struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
    
    typedef LinearAlgebra::RanderNorm<ScalarType,2> NormType;
    typedef LinearAlgebra::RanderNorm<ScalarType,1> NormType1;
    typedef NormType DistanceGuess;
};

struct TraitsAsymmetricQuadraticLag2 : TraitsBase<2> {
    typedef Lagrangian2Stencil<OffsetType,IndexDiff,Lagrangian2StencilPeriodicity::Simple> StencilType;
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
    Redeclare6Types(FromHFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,RecomputeType,Traits)
    Redeclare5Types(FromHFM,IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef)
    Redeclare6Types(FromTraits,NormType,NormType1,IndexType,StencilType,OffsetType,DistanceGuess)
    Redeclare1Type(FromSuperclass,OffsetVal3)
    Redeclare1Constant(FromHFM,Dimension)
    
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
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);
        pMetric = that->template GetField<MetricElementType>("metric",false);
    }
    virtual DistanceGuess GetGuess(IndexCRef index) const override;
private:
    std::forward_list<OffsetType> l;
    NormType GetNorm(IndexCRef index) const {
        const MetricElementType data = (*pMetric)(index);
        return NormType{data.first,data.second};}
};


typedef StencilQuadLinLag2<TraitsRanderLag2> StencilRanderLag2;
typedef StencilQuadLinLag2<TraitsAsymmetricQuadraticLag2> StencilAsymmetricQuadraticLag2;

#include "Lagrangian2Stencil.hxx"

#endif /* Lagrangian2Stencil_h */
