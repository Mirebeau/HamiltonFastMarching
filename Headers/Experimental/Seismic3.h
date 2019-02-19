//
//  Seismic3.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef Seismic3_h
#define Seismic3_h


#include "Base/Lagrangian3Stencil.h"
#include "JMM_CPPLibs/LinearAlgebra/SeismicNorm.h"
#include "Specializations/CommonTraits.h"

struct TraitsSeismic3 : TraitsBase<3> {
	using StencilType = Lagrangian3Stencil<OffsetType,ScalarType,DiscreteType>;
	using DomainType = PeriodicGrid<TraitsSeismic3>;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	using NormType = LinearAlgebra::SeismicNorm<ScalarType,3>;
	using DistanceGuess = NormType;
};


struct StencilSeismic3 : HamiltonFastMarching<TraitsSeismic3>::StencilDataType {
	typedef HamiltonFastMarching<TraitsSeismic3> HFM;
	typedef typename HFM::StencilDataType Superclass;
	
	Redeclare4Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType)
	Redeclare5Types(HFM,IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef)
	Redeclare4Types(HFM,RecomputeType,Traits,DomainType,IndexDiff)
	Redeclare5Types(TraitsSeismic3,NormType,IndexType,StencilType,OffsetType,DistanceGuess)
	Redeclare1Type(Superclass,OffsetVals)
	Redeclare1Constant(HFM,Dimension)
	
	// Specific to this model
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVals &) override;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;
	using MetricElementType = NormType::HookeTensorType;
	
	// Generic
	typedef typename Traits::template DataSource<MetricElementType> MetricType;
	std::unique_ptr<MetricType> pMetric;
	ParamDefault param;
	bool checkAcuteness = false;
	
	virtual void SetStencil(IndexCRef index, StencilType & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(IndexCRef index) const override {
		return GetNorm(index);}
private:
	std::forward_list<OffsetType> l;
	NormType GetNorm(IndexCRef index) const; // Includes rescaling by h
};

#include "Implementation/Seismic3.hpp"


#endif /* Seismic3_h */
