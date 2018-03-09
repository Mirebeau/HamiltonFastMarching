//
//  Riemannian45.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 15/02/2018.
//

#ifndef Riemannian45_h
#define Riemannian45_h

#include "Specializations/CommonTraits.h"
#include "LinearAlgebra/VoronoiReduction.h"

// ------------- 4D Riemannian metrics ------------
struct TraitsRiemann4 : TraitsBase<4> {
    typedef Difference<0> DifferenceType;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric = (Dimension*(Dimension+1))/2;
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsRiemann4::boundaryConditions) TraitsRiemann4::boundaryConditions;


struct StencilRiemann4 : HamiltonFastMarching<TraitsRiemann4>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann4> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    typedef VoronoiFirstReduction<ScalarType,Dimension> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],(*pDualMetric)(index));
        for(auto & diff : stencil.symmetric)
            diff.baseWeight/=square(param.gridScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");}
};

// ------------- 5D Riemannian metrics ------------
struct TraitsRiemann5 : TraitsBase<5> {
    typedef Difference<0> DifferenceType;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric = (Dimension*(Dimension+1))/2;
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsRiemann5::boundaryConditions) TraitsRiemann5::boundaryConditions;


struct StencilRiemann5 : HamiltonFastMarching<TraitsRiemann5>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann5> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    typedef VoronoiFirstReduction<ScalarType,Dimension> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],(*pDualMetric)(index));
        for(auto & diff : stencil.symmetric)
            diff.baseWeight/=square(param.gridScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");}
};

#endif /* Riemannian45_h */
