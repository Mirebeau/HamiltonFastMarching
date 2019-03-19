// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef RiemannLifted_h
#define RiemannLifted_h

#include "Specializations/CommonTraits.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"

template<Boundary cond>
struct TraitsRiemannLifted2 : TraitsBase<3> {
	using DifferenceType = EulerianDifference<OffsetType, ScalarType, 0>;
	using StencilType = EulerianStencil<DifferenceType, 3+1>; // symmetric
	
	typedef const std::array<Boundary,Dimension> BoundCondType;
	constexpr static BoundCondType boundaryConditions = {{Boundary::Closed, Boundary::Closed, cond}};
	using DomainType = PeriodicGrid<TraitsRiemannLifted2>;
	
    // Stencils actually depend on all coordinates. This is to get proper domain parametrization.
    static const DiscreteType nStencilDependencies=1;
    typedef const std::array<DiscreteType, nStencilDependencies> StencilDepType;
    constexpr static StencilDepType stencilDependencies = {{2}};
};
template<Boundary cond> constexpr typename TraitsRiemannLifted2<cond>::BoundCondType TraitsRiemannLifted2<cond>::boundaryConditions;
template<Boundary cond> constexpr typename TraitsRiemannLifted2<cond>::StencilDepType TraitsRiemannLifted2<cond>::stencilDependencies;


template<Boundary cond>
struct StencilRiemannLifted2 : HamiltonFastMarching<TraitsRiemannLifted2<cond> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemannLifted2<cond> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare7Types(Superclass,Traits,IndexType,StencilType,ParamInterface,HFMI,ScalarType,DiscreteType);
    typedef typename Traits::template BasisReduction<2> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType Sym;
    typedef LinearAlgebra::VectorPair<Sym,ScalarType> MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric, pMetric;
    
    typename HFM::ParamDefault param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        MetricElementType met;
        if(pDualMetric) met = (*pDualMetric)(index);
        else {met = (*pMetric)(index); met.first = met.first.Inverse(); met.second=1/met.second;}
		
		auto & symmetric = stencil.symmetric[0];
        Voronoi1Mat<ReductionType>(&symmetric[0], met.first/square(param.gridScale));
        symmetric[3].offset = {0,0,1};
        symmetric[3].baseWeight = met.second/square(param.dependScale);
    };
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        auto & io = that->io;
        const ScalarType bundleScale =
        io.template Get<ScalarType>("bundleScale",cond==Boundary::Periodic ? 2.*Traits::mathPi/this->dims[2] : 1.);
        param.Setup(that->io,bundleScale);
        if(io.HasField("metric")) pMetric = that->template GetField<MetricElementType>("metric");
        else pDualMetric=that->template GetField<MetricElementType>("dualMetric");
    }
};

struct TraitsRiemannLifted3 : TraitsBase<4> {
	using DifferenceType = EulerianDifference<OffsetType,ScalarType,0>;
	using StencilType = EulerianStencil<DifferenceType,6+1>;

	constexpr static std::array<Boundary,Dimension> boundaryConditions = {{Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed}};
	using DomainType = PeriodicGrid<TraitsRiemannLifted3>;
	
    // Stencils actually depend on all coordinates. This is to get proper domain parametrization.
    static const DiscreteType nStencilDependencies=1;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{3}};
};
constexpr decltype(TraitsRiemannLifted3::boundaryConditions) TraitsRiemannLifted3::boundaryConditions;
constexpr decltype(TraitsRiemannLifted3::stencilDependencies) TraitsRiemannLifted3::stencilDependencies;


struct StencilRiemannLifted3 : HamiltonFastMarching<TraitsRiemannLifted3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemannLifted3>::StencilDataType Superclass;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::SymmetricMatrixType Sym;
    typedef LinearAlgebra::VectorPair<Sym,ScalarType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    typename HFM::ParamDefault param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
		auto & symmetric = stencil.symmetric[0];
        Voronoi1Mat<ReductionType>(&symmetric[0], met.first/square(param.gridScale));
        symmetric[6].offset = {0,0,0,1};
        symmetric[6].baseWeight = met.second/square(param.dependScale);
    };
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        const ScalarType bundleScale = that->io.Get<ScalarType>("bundleScale",1.);
        param.Setup(that->io,bundleScale);
        pDualMetric=that->template GetField<MetricElementType>("dualMetric");}
};

#endif /* RiemannLifted_h */
