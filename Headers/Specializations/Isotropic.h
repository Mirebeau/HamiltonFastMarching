// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Isotropic_h
#define Isotropic_h

#include "CommonTraits.h"
#include "JMM_CPPLibs/LinearAlgebra/DiagonalNorm.h"

// ------------- Isotropic metrics ------------
template<size_t VDimension>
struct TraitsIsotropic : TraitsBase<VDimension> {
    using Superclass = TraitsBase<VDimension>;
    Redeclare3Types(Superclass,DiscreteType,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    using DifferenceType = EulerianDifference<OffsetType,ScalarType,1>;
    using StencilType = EulerianStencil<DifferenceType,Dimension>;
    using DomainType = PeriodicGrid<TraitsIsotropic>;
	using DistanceGuess = LinearAlgebra::DiagonalNorm<ScalarType,Dimension>;
};

template<size_t VDimension>
struct StencilIsotropic final :
HamiltonFastMarching<TraitsIsotropic<VDimension> >::StencilDataType {
    using HFM = HamiltonFastMarching<TraitsIsotropic<VDimension> >;
    using Superclass = typename HFM::StencilDataType;
    Redeclare10Types(HFM,ParamDefault,IndexType,StencilType,ParamInterface,
					HFMI,DistanceGuess,ScalarType,IndexCRef,PointType,VectorType)
    Redeclare1Constant(HFM,Dimension)
	using ParamType = typename HFM::template _ParamDefault<2,void>; // Distinct scale on each axis
	ParamType param;

    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        auto & differences = stencil.symmetric[0];
        for(int i=0; i<Dimension; ++i){
            auto & diff = differences[i];
            diff.baseWeight=1./square(param.gridScales[i]);
            diff.offset.fill(0);
            diff.offset[i]=1;
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {Superclass::Setup(that); param.Setup(that);}
    virtual DistanceGuess GetGuess(const PointType & p) const override {
		return Rescale(MapWeightedSum<ScalarType>(
				*this->pMultSource,this->pFM->dom.Neighbors(p)));}
	virtual DistanceGuess GetGuess(const IndexType & index) const override {
		return Rescale((*this->pMultSource)(index));}
protected:
	DistanceGuess Rescale(ScalarType s) const {
		DistanceGuess diag;
		const PointType h = param.gridScales;
		for(int i=0; i<Dimension; ++i) diag.d[i] = h[i]/s;
		return diag;
	}
};

// ----------- Diagonal metrics ------------

template<size_t VDimension>
struct TraitsDiagonal : TraitsBase<VDimension> {
    using Superclass = TraitsBase<VDimension>;
    Redeclare3Types(Superclass,DiscreteType,OffsetType,ScalarType)
    Redeclare1Constant(Superclass,Dimension)

    using DifferenceType = EulerianDifference<OffsetType,ScalarType,Dimension>;
    using StencilType = EulerianStencil<DifferenceType,Dimension>;
    using DomainType = PeriodicGrid<TraitsDiagonal>;
	using DistanceGuess = LinearAlgebra::DiagonalNorm<ScalarType,Dimension>;
};

template<size_t VDimension>
struct StencilDiagonal final
: HamiltonFastMarching<TraitsDiagonal<VDimension> >::StencilDataType {
    using HFM = HamiltonFastMarching<TraitsDiagonal<VDimension> >;
    using Superclass = typename HFM::StencilDataType;
    Redeclare4Types(HFM,IndexType,StencilType,ParamInterface,HFMI)
    Redeclare5Types(HFM,DistanceGuess,ScalarType,IndexCRef,VectorType,PointType)
    Redeclare1Constant(HFM,Dimension)
	using ParamType = typename HFM::template _ParamDefault<2,void>; // Distinct scale on each axis
    ParamType param;
    
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
		return GuessFromSpeeds((*this->pMultSource)(index));}
	virtual DistanceGuess GetGuess(const PointType & p) const override {
		auto vecSpeeds = [this](const IndexType & index)->VectorType {
			return VectorType::FromOrigin((*this->pMultSource)(index));};
		return GuessFromSpeeds(PointType::FromOrigin(
			MapWeightedSum<VectorType>(vecSpeeds,this->pFM->dom.Neighbors(p))));}
protected:
	DistanceGuess GuessFromSpeeds(const PointType & s) const {
		DistanceGuess diag;
		const PointType h = param.gridScales;
		for(int i=0; i<Dimension; ++i) diag.d[i] = h[i]/s[i];
		return diag;
	}
};

#endif /* Isotropic_h */
