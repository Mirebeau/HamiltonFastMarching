//
//  HalfDisk.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 10/07/2017.
//
//

#ifndef HalfDisk_h
#define HalfDisk_h

#include "Base/HamiltonFastMarching.h"
#include "Specializations/CommonTraits.h"
#include "LinearAlgebra/VectorPairType.h"

/*******
 An asymmetric metric, of "half disk" type.
 Metric is of dimension d+1, and is the concatenation of
 - a vector, indicating the speed forward.
 - a scalar, for sideways speed ratio.
 
 Relaxation parameters are:
 - a scalar eps, for how badly is backwards motion penalized.
 - a scalar epsForward, for how accurately is forward speed implemented.
 As usual, small parameters yield large stencils.
 
 
********/

// ------- Two dimensional model --------

struct TraitsHalfDisk2 : TraitsBase<2> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=3, nForward=3;
};

constexpr decltype(TraitsHalfDisk2::boundaryConditions) TraitsHalfDisk2::boundaryConditions;

struct StencilHalfDisk2 : HamiltonFastMarching<TraitsHalfDisk2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsHalfDisk2>::StencilDataType Superclass;
    ParamType param;
    ScalarType eps = 0.2, epsForward = 0.3;
    typedef Traits::BasisReduction<2> ReductionType;
    typedef LinearAlgebra::VectorPair<VectorType, ScalarType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;

    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        const VectorType & v = met.first/param.gridScale;
        const ScalarType & r = met.second;
        Voronoi1Vec<ReductionType>(&stencil.forward[0],
                                   v * sqrt(std::max(0.,1.-square(eps*r))),
                                   epsForward);
        typedef ReductionType::SymmetricMatrixType Sym;
        const Sym m =
        (Sym::RankOneTensor(v)*(square(eps)-1)
         +Sym::Identity()*v.SquaredNorm())
        *square(r);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m);
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        eps=that->io.Get<ScalarType>("eps",eps);
        epsForward=that->io.Get<ScalarType>("epsForward",eps*1.5);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");
    }
};


// ------- Three dimensional model ---------

struct TraitsHalfDisk3 : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=6, nForward=6;
};

constexpr decltype(TraitsHalfDisk3::boundaryConditions) TraitsHalfDisk3::boundaryConditions;

struct StencilHalfDisk3 : HamiltonFastMarching<TraitsHalfDisk3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsHalfDisk3>::StencilDataType Superclass;
    ParamType param;
    ScalarType eps = 0.2, epsForward = 0.3;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef LinearAlgebra::VectorPair<VectorType, ScalarType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        const VectorType & v = met.first/param.gridScale;
        const ScalarType & r = met.second;
        Voronoi1Vec<ReductionType>(&stencil.forward[0],
                                   v * sqrt(std::max(0.,1.-square(eps*r))),
                                   epsForward);
        typedef ReductionType::SymmetricMatrixType Sym;
        const Sym m =
        (Sym::RankOneTensor(v)*(square(eps)-1)
         +Sym::Identity()*v.SquaredNorm())
        *square(r);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m);
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        eps=that->io.Get<ScalarType>("eps",eps);
        epsForward=that->io.Get<ScalarType>("epsForward",eps*1.5);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");
    }
};

// ------- Three dimensional model with one additional radius dimension ---------

struct TraitsHalfDisk3p1 : TraitsBase<4> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=6+1, nForward=6;
};

constexpr decltype(TraitsHalfDisk3p1::boundaryConditions) TraitsHalfDisk3p1::boundaryConditions;

struct StencilHalfDisk3p1 : HamiltonFastMarching<TraitsHalfDisk3p1>::StencilDataType {
    typedef HamiltonFastMarching<TraitsHalfDisk3p1>::StencilDataType Superclass;
    ParamType param;
    ScalarType eps = 0.2, epsForward = 0.3;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::VectorType Vector3;
    /// Physical speed, bundle speed.
    typedef LinearAlgebra::Vector<ScalarType, 2> SpeedPair;
    typedef LinearAlgebra::VectorPair<Vector3, SpeedPair> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        const Vector3 & v = met.first/param.gridScale;
        const ScalarType & r = met.second[0];
        
        Voronoi1Vec<ReductionType>(&stencil.forward[0],
                                   v * sqrt(std::max(0.,1.-square(eps*r))),
                                   epsForward);
        typedef ReductionType::SymmetricMatrixType Sym;
        const Sym m =
        (Sym::RankOneTensor(v)*(square(eps)-1)
         +Sym::Identity()*v.SquaredNorm())
        *square(r);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m);
        
        for(int i=0; i<6; ++i){
            stencil.symmetric[i].offset[3]=0;
            stencil.forward[i].offset[3]=0;}
        
        stencil.symmetric[6].offset = OffsetType::Constant(0);
        stencil.symmetric[6].offset[3]=1;
        stencil.symmetric[6].baseWeight=square(met.second[1]/param.gridScale);
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        eps=that->io.Get<ScalarType>("eps",eps);
        epsForward=that->io.Get<ScalarType>("epsForward",eps*1.5);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");
    }
};


#endif /* HalfDisk_h */
