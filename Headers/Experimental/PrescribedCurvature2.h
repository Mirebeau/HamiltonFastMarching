// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef PrescribedCurvature2_h
#define PrescribedCurvature2_h

/**
 
 In this file, we implement minimal paths models which penalize the difference between the present value of path curvature, and some prescribed value. The present models generalize the ReedsShepp2, ReedsSheppForward2, Elastica2, and Dubins2 models defined in Curvature2.h.
 
 The speed function, prescribed curvature value, and the intensity of penalization of the curvature difference are set locally by the user. Due to this freedom, we cannot use a shared data structure for the stencils, hence the present models have a larger memory footprint and computation time than the original ones.
 
 The mathematical model for the metric reads, for a path parametrized at unit speed :
    F(x,x',x'') = C( xi (x''-kappa) )/speed,
 where speed (speed function), xi (curvature penalization intensity), kappa (prescribed curvature), are globally defined fields depending on x and x'.
 
 Time dependency, and Automatic differentiation, of the model parameters speed,xi,kappa, is presently not available. 
 //TODO : It is feasible for the speed function, but some architectural changes would be needed.
 
 */


#include "Specializations/Curvature2.h"


// --------------- R2S1NonShared Traits ---------------

struct TraitsR2S1NonShared : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Periodic}};
};

constexpr decltype(TraitsR2S1NonShared::boundaryConditions) TraitsR2S1NonShared::boundaryConditions;

// ----------- 2D Reeds-Shepp Extended model ---------
struct TraitsReedsSheppExt2 : TraitsR2S1NonShared {
    static const DiscreteType nSymmetric = 6+1;
};

struct StencilReedsSheppExt2
: HamiltonFastMarching<TraitsReedsSheppExt2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppExt2>::StencilDataType Superclass;
    typedef HamiltonFastMarching<TraitsR2S1>::StencilDataType::ParamType ParamType;
    ParamType param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pKappa;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index), xi=(*pXi)(index), kappa=(*pKappa)(index),
        gS=param.gridScale,tS=param.dependScale;
        const ScalarType theta = index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType v{c/gS,s/gS,kappa/tS};
        
        typedef Traits::BasisReduction<3> ReductionType;
        Voronoi1Vec<ReductionType>(&stencil.symmetric[0], speed*v, eps);
        
        stencil.symmetric[6].offset = OffsetType{0,0,1};
        stencil.symmetric[6].baseWeight = square(speed/(xi*tS));
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        pSpeed = that->GetField<ScalarType>("speed");
        pXi = that->GetField<ScalarType>("xi");
        pKappa = that->GetField<ScalarType>("kappa");
        param.Setup(that->io,2*mathPi/dims.back());}
};

// ----------- 2D Reeds-Shepp Forward Extended model ---------

struct TraitsReedsSheppForwardExt2 : TraitsR2S1NonShared {
    static const DiscreteType nSymmetric = 1, nForward=6;
};

struct StencilReedsSheppForwardExt2
: HamiltonFastMarching<TraitsReedsSheppForwardExt2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForwardExt2>::StencilDataType Superclass;
    typedef HamiltonFastMarching<TraitsR2S1>::StencilDataType::ParamType ParamType;
    ParamType param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pKappa;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index), xi=(*pXi)(index), kappa=(*pKappa)(index),
        gS=param.gridScale,tS=param.dependScale;
        const ScalarType theta = index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType v{c/gS,s/gS,kappa/tS};
        
        typedef Traits::BasisReduction<3> ReductionType;
        Voronoi1Vec<ReductionType>(&stencil.forward[0], speed*v, eps);
        
        stencil.symmetric[0].offset = OffsetType{0,0,1};
        stencil.symmetric[0].baseWeight = square(speed/(xi*tS));
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        pSpeed = that->GetField<ScalarType>("speed");
        pXi = that->GetField<ScalarType>("xi");
        pKappa = that->GetField<ScalarType>("kappa");
        param.Setup(that->io,2*mathPi/dims.back());}
};


// ---------- 2D Dubins Extended model -------------

struct TraitsDubinsExt2 : TraitsR2S1NonShared {
    static const DiscreteType nMax=2, nMaxForward=6;
};

struct StencilDubinsExt2
: HamiltonFastMarching<TraitsDubinsExt2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDubinsExt2>::StencilDataType Superclass;
    typedef HamiltonFastMarching<TraitsR2S1>::StencilDataType::ParamType ParamType;
    ParamType param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pKappa;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index), xi=(*pXi)(index), kappa=(*pKappa)(index),
        gS=param.gridScale,tS=param.dependScale;
        const ScalarType theta = index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType
        vL{c/gS,s/gS,(kappa+1./xi)/tS},
        vR{c/gS,s/gS,(kappa-1./xi)/tS};
        
        typedef Traits::BasisReduction<3> ReductionType;
        Voronoi1Vec<ReductionType>(&stencil.maxForward[0][0], speed*vL, eps);
        Voronoi1Vec<ReductionType>(&stencil.maxForward[1][0], speed*vR, eps);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        pSpeed = that->GetField<ScalarType>("speed");
        pXi = that->GetField<ScalarType>("xi");
        pKappa = that->GetField<ScalarType>("kappa");
        param.Setup(that->io,2*mathPi/dims.back());}
};


// ----------- 2D Euler-Mumford Elastica Extended model -----------

template<int nFejer>
struct TraitsElasticaExt2 : TraitsR2S1NonShared {
    static const DiscreteType nForward=6*nFejer;
};

template<int nFejer>
struct StencilElasticaExt2
: HamiltonFastMarching<TraitsElasticaExt2<nFejer> >::StencilDataType {
    typedef typename HamiltonFastMarching<TraitsElasticaExt2<nFejer> >::StencilDataType Superclass;
    Redeclare7Types(FromSuperclass,Traits,ScalarType,IndexCRef,VectorType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(FromTraits,mathPi)
    typedef HamiltonFastMarching<TraitsR2S1>::StencilDataType::ParamType ParamType;
    ParamType param;
    ScalarType eps=0.1;
    typedef typename Traits::template DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pKappa;

    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index), xi=(*pXi)(index), kappa=(*pKappa)(index),
        gS=param.gridScale,tS=param.dependScale;
        const ScalarType theta = index[2]*tS;
        const ScalarType cT = cos(theta), sT=sin(theta);
        
        for(int l=0; l<nFejer; ++l){
            const ScalarType phi = mathPi*(l+0.5)/nFejer;
            const ScalarType cP = cos(phi), sP=sin(phi);
            const VectorType v{sP*cT/gS,sP*sT/gS,(sP*kappa+cP/xi)/tS};
            typedef typename Traits::template BasisReduction<3> ReductionType;
            Voronoi1Vec<ReductionType>(&stencil.forward[6*l], v*speed, eps);
            for(int i=0; i<6; ++i) stencil.forward[6*l+i].baseWeight*=StencilElastica2<nFejer>::fejerWeights[l];
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        pSpeed = that->template GetField<ScalarType>("speed");
        pXi = that->template GetField<ScalarType>("xi");
        pKappa = that->template GetField<ScalarType>("kappa");
        param.Setup(that->io,2*mathPi/this->dims.back());}
};



#endif /* PrescribedCurvature2_h */
