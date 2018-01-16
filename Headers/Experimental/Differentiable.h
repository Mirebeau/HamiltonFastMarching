// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2018 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Differentiable_h
#define Differentiable_h

/* In this file, we define differentiable implementations of some common models.
 More precisely we replace symmetric upwind finite differences
    max(0,U(x)-U(x-e),U(x)-U(x+e))^2
 with two upwind finite differences
    max(0,U(x)-U(x-e))^2 + max(0,U(x)-U(x+e))^2
 
 This new implementation is expected to be slightly less accurate near the cut locus, and marginally slower.
 However the value function should be continuously differentiable w.r.t. the parameters (cost, boundary conditions),
 which may prove useful in optimization tasks.
 */

#include "Specializations/CommonTraits.h"
#include "Base/HFMInterface.h"


// ------------- 2D Isotropic metric, differentiable implementation ------------
struct TraitsIsotropicDiff2 : TraitsBase<2> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{}};
    static const DiscreteType nForward = 2*Dimension;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions = {{Boundary::Closed, Boundary::Closed}};
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsIsotropicDiff2::stencilDependencies) TraitsIsotropicDiff2::stencilDependencies;
constexpr const decltype(TraitsIsotropicDiff2::boundaryConditions)  TraitsIsotropicDiff2::boundaryConditions;

struct StencilIsotropicDiff2 : HamiltonFastMarching<TraitsIsotropicDiff2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropicDiff2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        int n=0;
        for(int i=0; i<Dimension; ++i){
            for(int s=-1; s<=1; s+=2){
                auto & diff = stencil.forward[n]; ++n;
                diff.baseWeight=1./square(param.gridScale);
                diff.offset.fill(0);
                diff.offset[i]=s;
            }
        }
        assert(n==stencil.forward.size());
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};

#endif /* Differentiable_h */
