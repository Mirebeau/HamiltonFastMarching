// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#include "CommonTraits.h"

// ------------- 2D - 3D Riemannian metrics ------------
// Bounded box in all dimensions, except possibly the last one.
template<size_t VDimension, Boundary lastBoundary>
struct TraitsRiemann : TraitsBase<VDimension> {
    using Superclass = TraitsBase<VDimension>;
    Redeclare2Types(Superclass,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    using DifferenceType = EulerianDifference<OffsetType,ScalarType,0>;
    using StencilType = EulerianStencil<DifferenceType,(Dimension*(Dimension+1))/2>;

	constexpr static const Boundary_ClosedButLast<Dimension, lastBoundary> boundaryConditions{};
    using DomainType = PeriodicGrid<TraitsRiemann>;
};

template<size_t VDimension, Boundary lastBoundary = Boundary::Closed>
struct StencilRiemann final
: HamiltonFastMarching<TraitsRiemann<VDimension,lastBoundary> >::StencilDataType {
    using HFM = HamiltonFastMarching<TraitsRiemann<VDimension,lastBoundary> >;
    using Superclass = typename HFM::StencilDataType ;
    Redeclare10Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI,Traits,
					ScalarType,IndexCRef,DistanceGuess,PointType)
    Redeclare1Constant(HFM,Dimension)
	using ParamType = typename HFM::template _ParamDefault<2,void>; // Distinct scale on each axis
    ParamType param;

    using ReductionType = typename Traits::template BasisReduction<Dimension>;
    using SymmetricMatrixType = typename ReductionType::SymmetricMatrixType;
    using MetricElementType = SymmetricMatrixType;
    using MetricType = typename Traits::template DataSource<MetricElementType>;
    std::unique_ptr<MetricType> pDualMetric;
    std::unique_ptr<MetricType> pMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
		const SymmetricMatrixType diff =
		Rescale<true>(pDualMetric ? (*pDualMetric)(index) : (*pMetric)(index).Inverse() );
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0][0], diff);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
		for(int i=0; i<Dimension; ++i) for(int j=0;j<=i;++j) {
			auto & h = param.gridScales;
			gridScales(i,j)=h[i]*h[j];
			invGridScales(i,j)=1/gridScales(i,j);
		}
        if(that->io.HasField("dualMetric")) pDualMetric = that->template GetField<MetricElementType>("dualMetric");
        else pMetric = that->template GetField<MetricElementType>("metric");
    }
    virtual DistanceGuess GetGuess(const PointType & p) const override {
		return Rescale<false>(pDualMetric ?
		MapWeightedSum<SymmetricMatrixType>(*pDualMetric,this->pFM->dom.Neighbors(p)).Inverse() :
		MapWeightedSum<SymmetricMatrixType>(*pMetric,this->pFM->dom.Neighbors(p)) );
	}
	virtual DistanceGuess GetGuess(const IndexType & index) const override {
		return Rescale<false>(pDualMetric ? (*pDualMetric)(index).Inverse() : (*pMetric)(index));
	}
protected:
	SymmetricMatrixType gridScales, invGridScales; // Stores h[i]*h[j] where h is the gridscale
	template<bool dual> SymmetricMatrixType Rescale(SymmetricMatrixType m) const {
		for(int i=0;i<Dimension;++i) for(int j=0;j<=i;++j) m(i,j)*=(dual ? invGridScales(i,j) : gridScales(i,j));
		return m;}
};

