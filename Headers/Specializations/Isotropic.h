// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Isotropic_h
#define Isotropic_h

#include "CommonTraits.h"
#include "Base/HFMInterface.h"


// ------------- Isotropic metrics ------------
template<size_t VDimension>
struct TraitsIsotropic : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare1Type(FromSuperclass,DiscreteType)
    Redeclare1Constant(FromSuperclass,Dimension)
    
    typedef typename Superclass::template Difference<1> DifferenceType;
//    static const DiscreteType nStencilDependencies=0;
//    typedef std::array<DiscreteType, nStencilDependencies> StencilDepType;
//    constexpr static const StencilDepType stencilDependencies = {{}};
    
    static const DiscreteType nSymmetric = Dimension;
//    constexpr static const Boundary_AllClosed boundaryConditions{};
};
// Linker wants the following line for some obscure reason.
//template<size_t VD> constexpr const typename TraitsIsotropic<VD>::StencilDepType TraitsIsotropic<VD>::stencilDependencies;
//template<size_t VD> constexpr const Boundary_AllClosed TraitsIsotropic<VD>::boundaryConditions;

template<size_t VDimension>
struct StencilIsotropic : HamiltonFastMarching<TraitsIsotropic<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropic<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(FromHFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(FromHFM,Dimension)
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
    Redeclare1Type(FromSuperclass,DiscreteType)
    Redeclare1Constant(FromSuperclass,Dimension)

    typedef typename Superclass::template Difference<Dimension> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    typedef std::array<DiscreteType, nStencilDependencies> StencilDepType;
    constexpr static const StencilDepType stencilDependencies = {{}};
    
    static const DiscreteType nSymmetric = Dimension;
    constexpr static const Boundary_AllClosed boundaryConditions{};
};
template<size_t VD> constexpr const typename TraitsDiagonal<VD>::StencilDepType TraitsDiagonal<VD>::stencilDependencies;
template<size_t VD> constexpr const Boundary_AllClosed TraitsDiagonal<VD>::boundaryConditions;

template<size_t VDimension>
struct StencilDiagonal : HamiltonFastMarching<TraitsDiagonal<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsDiagonal<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare4Types(FromHFM,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(FromHFM,Dimension)
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


// --------- Previous non-templated instantiations ---------

/*
struct TraitsIsotropic2 : TraitsBase<2> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{}};
    static const DiscreteType nSymmetric = Dimension;
    constexpr static const Boundary_AllClosed boundaryConditions{};
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsIsotropic2::stencilDependencies) TraitsIsotropic2::stencilDependencies;
constexpr const decltype(TraitsIsotropic2::boundaryConditions) TraitsIsotropic2::boundaryConditions;

struct StencilIsotropic2 : HamiltonFastMarching<TraitsIsotropic2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropic2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    
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
*/


// --- 2D diagonal metrics ---
/*
struct TraitsDiagonal2 : TraitsBase<2> {
    typedef Difference<2> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{}};
    static const DiscreteType nSymmetric = Dimension;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions = {{Boundary::Closed, Boundary::Closed}};
};
constexpr const decltype(TraitsDiagonal2::boundaryConditions) TraitsDiagonal2::boundaryConditions;
constexpr const decltype(TraitsDiagonal2::stencilDependencies) TraitsDiagonal2::stencilDependencies;

struct StencilDiagonal2 : HamiltonFastMarching<TraitsDiagonal2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDiagonal2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::_ParamDefault<2> param;
    
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
};*/

// -- 3D isotropic metrics ---
// Cannot template over dimension, because 'boundaryConditions' array is otherwise impossible to initialize

/*
struct TraitsIsotropic3 : TraitsBase<3> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{}};
    static const DiscreteType nSymmetric = Dimension;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions = {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
};
constexpr const decltype(TraitsIsotropic3::stencilDependencies) TraitsIsotropic3::stencilDependencies;
constexpr const decltype(TraitsIsotropic3::boundaryConditions) TraitsIsotropic3::boundaryConditions;

struct StencilIsotropic3 : HamiltonFastMarching<TraitsIsotropic3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropic3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    
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
*/
// --- 3D diagonal metrics ---
/*
struct TraitsDiagonal3 : TraitsBase<3> {
    typedef Difference<3> DifferenceType;
    static const DiscreteType nStencilDependencies=0;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{}};
    static const DiscreteType nSymmetric = Dimension;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions = {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
};
constexpr const decltype(TraitsDiagonal3::boundaryConditions) TraitsDiagonal3::boundaryConditions;
constexpr const decltype(TraitsDiagonal3::stencilDependencies) TraitsDiagonal3::stencilDependencies;

struct StencilDiagonal3 : HamiltonFastMarching<TraitsDiagonal3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDiagonal3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::_ParamDefault<2> param;
    
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
*/
#endif /* Isotropic_h */
