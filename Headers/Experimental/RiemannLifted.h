// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef RiemannLifted_h
#define RiemannLifted_h

#include "Specializations/CommonTraits.h"
#include "Base/HamiltonFastMarching.h"
#include "LinearAlgebra/VectorPairType.h"

//template<typename HFM> struct ParamLifted : typename

template<Boundary cond>
struct TraitsRiemannLifted2 : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    typedef const std::array<Boundary,Dimension> BoundCondType;
    constexpr static BoundCondType boundaryConditions = {{Boundary::Closed, Boundary::Closed, cond}};
    static const DiscreteType nSymmetric=3+1;
};
template<Boundary cond> constexpr typename TraitsRiemannLifted2<cond>::BoundCondType TraitsRiemannLifted2<cond>::boundaryConditions;


template<Boundary cond>
struct StencilRiemannLifted2 : HamiltonFastMarching<TraitsRiemannLifted2<cond> >::StencilDataType {
    typedef typename HamiltonFastMarching<TraitsRiemannLifted2<cond> >::StencilDataType Superclass;
    Redeclare7Types(FromSuperclass,Traits,IndexType,StencilType,ParamInterface,HFMI,ScalarType,DiscreteType);
    typedef typename Traits::template BasisReduction<2> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType Sym;
    typedef LinearAlgebra::VectorPair<Sym,ScalarType> MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    // Traits used to get a specific scale for the last coordinate
    struct DummyTraits : Traits {
        typedef typename Traits::template Difference<1> DifferenceType;
        static const DiscreteType nStencilDependencies = 1;
        typedef const std::array<DiscreteType, nStencilDependencies> StencilDepType;
        constexpr static StencilDepType stencilDependencies = {{2}};};
    typedef typename HamiltonFastMarching<DummyTraits>::StencilDataType::ParamType ParamType;
    ParamType param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], met.first/square(param.gridScale));
        stencil.symmetric[3].offset = {0,0,1};
        stencil.symmetric[3].baseWeight = met.second/square(param.dependScale);
    };
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        const ScalarType bundleScale =
        that->io.template Get<ScalarType>("bundleScale",cond==Boundary::Periodic ? 2.*Traits::mathPi/this->dims[2] : 1.);
        param.Setup(that->io,bundleScale);
        pDualMetric=that->template GetField<MetricElementType>("dualMetric");        
    }
};
template<Boundary cond> constexpr typename StencilRiemannLifted2<cond>::DummyTraits::StencilDepType
StencilRiemannLifted2<cond>::DummyTraits::stencilDependencies;


struct TraitsRiemannLifted3 : TraitsBase<4> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary,Dimension> boundaryConditions = {{Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=6+1;
};
constexpr decltype(TraitsRiemannLifted3::boundaryConditions) TraitsRiemannLifted3::boundaryConditions;


struct StencilRiemannLifted3 : HamiltonFastMarching<TraitsRiemannLifted3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemannLifted3>::StencilDataType Superclass;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::SymmetricMatrixType Sym;
    typedef LinearAlgebra::VectorPair<Sym,ScalarType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    // Traits used to get a specific scale for the last coordinate
    struct DummyTraits : Traits {
        typedef Difference<1> DifferenceType; static const DiscreteType nStencilDependencies = 1;
        constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{3}};};
    typedef HamiltonFastMarching<DummyTraits>::StencilDataType::ParamType ParamType;
    ParamType param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], met.first/square(param.gridScale));
        stencil.symmetric[6].offset = {0,0,0,1};
        stencil.symmetric[6].baseWeight = met.second/square(param.dependScale);
    };
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        const ScalarType bundleScale = that->io.Get<ScalarType>("bundleScale",1.);
        param.Setup(that->io,bundleScale);
        pDualMetric=that->template GetField<MetricElementType>("dualMetric");}
};
constexpr decltype(StencilRiemannLifted3::DummyTraits::stencilDependencies) StencilRiemannLifted3::DummyTraits::stencilDependencies;


#endif /* RiemannLifted_h */
