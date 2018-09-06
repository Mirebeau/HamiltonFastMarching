// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef IsotropicSpecializations_h
#define IsotropicSpecializations_h

#include "CommonTraits.h"

// -- Isotropic metrics on 2D and 3D domains with various boundary conditions (closed, periodic, klein bottle) --
// Used mainly for testing these boundary conditions.

template<Boundary cond>
struct TraitsIsotropicBox2 : TraitsBase<2> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    typedef const std::array<DiscreteType, nStencilDependencies> StencilDepType;
    typedef const std::array<Boundary, Dimension> BoundCondType;
    
    constexpr static StencilDepType stencilDependencies = {};
    constexpr static BoundCondType boundaryConditions =
    {{cond, cond==Boundary::Sphere2_0 ? Boundary::Sphere2_1 : cond}};
    static const DiscreteType nSymmetric=Dimension ;
};
template<Boundary cond>  constexpr typename TraitsIsotropicBox2<cond>::StencilDepType TraitsIsotropicBox2<cond>::stencilDependencies;
template<Boundary cond>  constexpr typename TraitsIsotropicBox2<cond>::BoundCondType TraitsIsotropicBox2<cond>::boundaryConditions;


template<Boundary cond>
struct StencilIsotropicBox2
: HamiltonFastMarching<TraitsIsotropicBox2<cond> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropicBox2<cond> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare4Types(Superclass,IndexType,StencilType,ParamInterface,HFMI)
    typedef typename HFM::ParamDefault ParamType;
    ParamType param;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        int i=0;
        for(auto & diff : stencil.symmetric){
            diff.baseWeight=1./square(param.gridScale);
            diff.offset.fill(0);
            diff.offset[i++]=1;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};

// ------- 2D Sphere ------
struct TraitsSphere2 : TraitsBase<2> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies = 1;
    constexpr static std::array<DiscreteType, nStencilDependencies>
        stencilDependencies = {{0}};
    constexpr static std::array<Boundary,Dimension> boundaryConditions =
        {{Boundary::Sphere2_0, Boundary::Sphere2_1}};
    static const DiscreteType nSymmetric=Dimension;
};
constexpr decltype(TraitsSphere2::stencilDependencies) TraitsSphere2::stencilDependencies;
constexpr decltype(TraitsSphere2::boundaryConditions) TraitsSphere2::boundaryConditions;

struct StencilSphere2 : HamiltonFastMarching<TraitsSphere2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsSphere2> HFM;
    typedef HFM::StencilDataType Superclass;
    typedef typename HFM::ParamDefault ParamType;
    ParamType param;
    bool projective = false;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        auto & diff0 = stencil.symmetric[0];
        diff0.baseWeight=1./square(param.gridScale);
        diff0.offset={1,0};
        
        const ScalarType theta =
            (index[0]+0.5)*(projective ? 0.5*mathPi : mathPi)/dims[0];
        auto & diff1 = stencil.symmetric[1];
        diff1.baseWeight=1./square(sin(theta)*param.gridScale);
        diff1.offset={0,1};
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that);
    projective = (bool)that->io.template Get<ScalarType>("projective",0.);
        param.origin.back() = -param.gridScale/2;
        if(( projective && dims[1]!=4*dims[0]) ||
           (!projective && dims[1]!=2*dims[0]))
            Msg() << "Sphere will not be round.\n";
    }
};


// ---------- 3D box domain with various boundary conditions ----------


template<Boundary cond>
struct TraitsIsotropicBox3 : TraitsBase<3> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    typedef const std::array<DiscreteType, nStencilDependencies> StencilDepType;
    typedef const std::array<Boundary, Dimension> BoundCondType;
    
    constexpr static StencilDepType stencilDependencies = {};
    constexpr static BoundCondType boundaryConditions =
    {{cond, cond, cond==Boundary::Sphere2_0 ? Boundary::Sphere2_1 : cond}};
    static const DiscreteType nSymmetric=Dimension ;
};
template<Boundary cond>  constexpr typename TraitsIsotropicBox3<cond>::StencilDepType TraitsIsotropicBox3<cond>::stencilDependencies;
template<Boundary cond>  constexpr typename TraitsIsotropicBox3<cond>::BoundCondType TraitsIsotropicBox3<cond>::boundaryConditions;


template<Boundary cond>
struct StencilIsotropicBox3
: HamiltonFastMarching<TraitsIsotropicBox3<cond> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropicBox3<cond> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare4Types(Superclass,IndexType,StencilType,ParamInterface,HFMI)
    typedef typename HFM::ParamDefault ParamType;
    ParamType param;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        int i=0;
        for(auto & diff : stencil.symmetric){
            diff.baseWeight=1./square(param.gridScale);
            diff.offset.fill(0);
            diff.offset[i++]=1;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};

#endif /* IsotropicSpecializations_h */
