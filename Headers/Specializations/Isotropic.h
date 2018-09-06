// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Isotropic_h
#define Isotropic_h

#include "CommonTraits.h"

// ------------- Isotropic metrics ------------
template<size_t VDimension>
struct TraitsIsotropic : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare1Type(Superclass,DiscreteType)
    Redeclare1Constant(Superclass,Dimension)
    typedef typename Superclass::template Difference<1> DifferenceType;
    static const DiscreteType nSymmetric = Dimension;
};

template<size_t VDimension>
struct StencilIsotropic : HamiltonFastMarching<TraitsIsotropic<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropic<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        for(int i=0; i<Dimension; ++i){
            auto & diff = stencil.symmetric[i];
            diff.baseWeight=1./square(param.gridScale);
            diff.offset.fill(0);
            diff.offset[i]=1;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};

// ----------- Diagonal metrics ------------

template<size_t VDimension>
struct TraitsDiagonal : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare1Type(Superclass,DiscreteType)
    Redeclare1Constant(Superclass,Dimension)

    typedef typename Superclass::template Difference<Dimension> DifferenceType;
    static const DiscreteType nSymmetric = Dimension;
};

template<size_t VDimension>
struct StencilDiagonal : HamiltonFastMarching<TraitsDiagonal<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsDiagonal<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare4Types(HFM,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(HFM,Dimension)
    typename HFM::template _ParamDefault<2> param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        for(int i=0; i<Dimension; ++i){
            auto & diff = stencil.symmetric[i];
            diff.baseWeight=1./square(param.gridScales[i]);
            diff.offset.fill(0);
            diff.offset[i]=1;
            diff.multIndex = i;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
};

#endif /* Isotropic_h */
