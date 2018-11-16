// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#include "CommonTraits.h"

// ------------- 2D - 3D Riemannian metrics ------------
template<size_t VDimension>
struct TraitsRiemann : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare2Types(Superclass,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    typedef EulerianDifference<OffsetType,ScalarType,0> DifferenceType;
    typedef EulerianStencil<DifferenceType,(Dimension*(Dimension+1))/2> StencilType;

    typedef PeriodicGrid<TraitsRiemann> DomainType;
};

template<size_t VDimension>
struct StencilRiemann : HamiltonFastMarching<TraitsRiemann<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare6Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI,Traits)
    Redeclare3Types(HFM,ScalarType,IndexCRef,DistanceGuess)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;

    typedef typename Traits::template BasisReduction<Dimension> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    std::unique_ptr<MetricType> pMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0][0],
                                   pDualMetric ? (*pDualMetric)(index) : (*pMetric)(index).Inverse() );
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff : stencil.symmetric[0]) diff.baseWeight*=hm2;
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        if(that->io.HasField("dualMetric")) pDualMetric = that->template GetField<MetricElementType>("dualMetric");
        else pMetric = that->template GetField<MetricElementType>("metric");
    }
    virtual DistanceGuess GetGuess(IndexCRef index) const override {
        const SymmetricMatrixType m = pDualMetric ? (*pDualMetric)(index).Inverse() : (*pMetric)(index);
        return m*square(param.gridScale);}
};

