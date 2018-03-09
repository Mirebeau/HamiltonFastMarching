// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#include "CommonTraits.h"
#include "Base/HFMInterface.h"


// ------------- 2D Riemannian metrics ------------
struct TraitsRiemann2 : TraitsBase<2> {
    typedef Difference<0> DifferenceType;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric = 3;
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsRiemann2::boundaryConditions) TraitsRiemann2::boundaryConditions;


struct StencilRiemann2 : HamiltonFastMarching<TraitsRiemann2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    typedef Traits::BasisReduction<2> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    std::unique_ptr<MetricType> pMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],
                                   pDualMetric ? (*pDualMetric)(index) : (*pMetric)(index).Inverse() );
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff : stencil.symmetric) diff.baseWeight*=hm2;
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        if(that->io.HasField("dualMetric")) pDualMetric = that->GetField<MetricElementType>("dualMetric");
        else pMetric = that->GetField<MetricElementType>("metric");
    }
};

// --------------- 3D Riemannian metrics ------------
struct TraitsRiemann3 : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric = 6;
};
constexpr decltype(TraitsRiemann3::boundaryConditions) TraitsRiemann3::boundaryConditions;

struct StencilRiemann3 : HamiltonFastMarching<TraitsRiemann3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRiemann3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef SymmetricMatrixType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],
                                   pDualMetric ? (*pDualMetric)(index) : (*pMetric)(index).Inverse() );
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff : stencil.symmetric) diff.baseWeight*=hm2;
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        if(that->io.HasField("dualMetric")) pDualMetric = that->GetField<MetricElementType>("dualMetric");
        else pMetric = that->GetField<MetricElementType>("metric");
    }};
