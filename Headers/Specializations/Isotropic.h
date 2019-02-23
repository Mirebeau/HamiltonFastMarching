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
    Redeclare3Types(Superclass,DiscreteType,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    typedef EulerianDifference<OffsetType,ScalarType,1> DifferenceType;
    typedef EulerianStencil<DifferenceType,Dimension> StencilType;
    typedef PeriodicGrid<TraitsIsotropic> DomainType;
};

template<size_t VDimension>
struct StencilIsotropic final :
HamiltonFastMarching<TraitsIsotropic<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsIsotropic<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare4Types(HFM,DistanceGuess,ScalarType,IndexCRef,PointType)
    Redeclare1Constant(HFM,Dimension)
    ParamDefault param;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        auto & differences = stencil.symmetric[0];
        for(int i=0; i<Dimension; ++i){
            auto & diff = differences[i];
            diff.baseWeight=1./square(param.gridScale);
            diff.offset.fill(0);
            diff.offset[i]=1;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
    virtual DistanceGuess GetGuess(const PointType & p) const override final {
		const ScalarType s = MapWeightedSum<ScalarType>(*this->pMultSource,this->pFM->dom.Neighbors(p));
        const ScalarType h = param.gridScale;
        return DistanceGuess::Identity() * square(h/s);}
};

// ----------- Diagonal metrics ------------

template<size_t VDimension>
struct TraitsDiagonal : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare3Types(Superclass,DiscreteType,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    typedef EulerianDifference<OffsetType,ScalarType,Dimension> DifferenceType;
    typedef EulerianStencil<DifferenceType,Dimension> StencilType;
    
    typedef PeriodicGrid<TraitsDiagonal> DomainType;
};

template<size_t VDimension>
struct StencilDiagonal final
: HamiltonFastMarching<TraitsDiagonal<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsDiagonal<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare4Types(HFM,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare5Types(HFM,DistanceGuess,ScalarType,IndexCRef,VectorType,PointType)
    Redeclare1Constant(HFM,Dimension)
    typename HFM::template _ParamDefault<2,void> param;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        for(int i=0; i<Dimension; ++i){
            auto & diff = stencil.symmetric[0][i];
            diff.baseWeight=1./square(param.gridScales[i]);
            diff.offset.fill(0);
            diff.offset[i]=1;
            diff.multIndex = i;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
    virtual DistanceGuess GetGuess(IndexCRef index) const override {
        const PointType h = param.gridScales;
        const PointType s=(*this->pMultSource)(index);
        VectorType diag; for(int i=0; i<Dimension; ++i) diag[i] = square(h[i]/s[i]);
        return DistanceGuess::Diagonal(diag);}
};

#endif /* Isotropic_h */
