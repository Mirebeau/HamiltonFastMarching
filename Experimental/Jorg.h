// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Jorg_h
#define Jorg_h

#include "Specializations/CommonTraits.h"
#include "Base/HamiltonFastMarching.h"
#include "LinearAlgebra/VectorPairType.h"

struct TraitsAdaptiveJorg : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Periodic}};
    static const DiscreteType nSymmetric=3+1;
};
constexpr decltype(TraitsAdaptiveJorg::boundaryConditions) TraitsAdaptiveJorg::boundaryConditions;


struct StencilAdaptiveJorg : HamiltonFastMarching<TraitsAdaptiveJorg>::StencilDataType {
    typedef HamiltonFastMarching<TraitsAdaptiveJorg>::StencilDataType Superclass;
    ParamType param;
    typedef Traits::BasisReduction<2> ReductionType;
    typedef ReductionType::SymmetricMatrixType Sym;
    typedef LinearAlgebra::VectorPair<Sym,ScalarType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], met.first);
        stencil.symmetric[3].offset = {0,0,1};
        stencil.symmetric[3].baseWeight = met.second;
        for(auto & diff : stencil.symmetric)
            diff.baseWeight/=square(param.gridScale);
    };
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);
        pDualMetric=that->template GetField<MetricElementType>("dualMetric");}
};

#endif /* Jorg_h */
