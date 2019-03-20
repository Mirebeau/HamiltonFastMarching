// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#include "CommonTraits.h"

// ------------- 2D - 3D Riemannian metrics ------------
// Bounded box in all dimensions, except possibly the last one.
template<size_t VDimension, Boundary lastBoundary>
struct TraitsRiemann : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare2Types(Superclass,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    typedef EulerianDifference<OffsetType,ScalarType,0> DifferenceType;
    typedef EulerianStencil<DifferenceType,(Dimension*(Dimension+1))/2> StencilType;

	constexpr static const Boundary_ClosedButLast<Dimension, lastBoundary> boundaryConditions{};
    typedef PeriodicGrid<TraitsRiemann> DomainType;
};

template<size_t VDimension, Boundary lastBoundary = Boundary::Closed>
struct StencilRiemann final
: HamiltonFastMarching<TraitsRiemann<VDimension,lastBoundary> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann<VDimension,lastBoundary> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare10Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI,Traits,
					ScalarType,IndexCRef,DistanceGuess,PointType)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;

    typedef typename Traits::template BasisReduction<Dimension> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    std::unique_ptr<MetricType> pMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
		const SymmetricMatrixType diff = (1./square(param.gridScale)) *
		(pDualMetric ? (*pDualMetric)(index) : (*pMetric)(index).Inverse() );
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0][0], diff);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        if(that->io.HasField("dualMetric")) pDualMetric = that->template GetField<MetricElementType>("dualMetric");
        else pMetric = that->template GetField<MetricElementType>("metric");
    }
    virtual DistanceGuess GetGuess(const PointType & p) const override {
		return square(param.gridScale) * (pDualMetric ?
		MapWeightedSum<SymmetricMatrixType>(*pDualMetric,this->pFM->dom.Neighbors(p)).Inverse() :
		MapWeightedSum<SymmetricMatrixType>(*pMetric,this->pFM->dom.Neighbors(p)) );
	}
	virtual DistanceGuess GetGuess(const IndexType & index) const override {
		return square(param.gridScale) *
		(pDualMetric ? (*pDualMetric)(index).Inverse() : (*pMetric)(index));
	}
};

