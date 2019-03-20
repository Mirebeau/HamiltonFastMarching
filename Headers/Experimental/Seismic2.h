//
//  Seismic2.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 30/01/2019.
//

#ifndef Seismic2_h
#define Seismic2_h

#include "Base/Lagrangian2Stencil.h"
#include "JMM_CPPLibs/LinearAlgebra/SeismicNorm.h"
#include "Specializations/CommonTraits.h"

struct TraitsSeismic2 : TraitsBase<2> {
	typedef Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType> StencilType;
	typedef PeriodicGrid<TraitsSeismic2> DomainType;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	typedef LinearAlgebra::SeismicNorm<ScalarType,2> NormType;
	typedef NormType DistanceGuess;
};


struct StencilSeismic2 final
: HamiltonFastMarching<TraitsSeismic2>::StencilDataType {
	typedef HamiltonFastMarching<TraitsSeismic2> HFM;
	typedef typename HFM::StencilDataType Superclass;

	Redeclare13Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,
					IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,RecomputeType,
					Traits,DomainType,IndexDiff)
	Redeclare5Types(TraitsSeismic2,NormType,IndexType,StencilType,OffsetType,DistanceGuess)
	Redeclare1Type(Superclass,OffsetVal3)
	Redeclare1Constant(HFM,Dimension)
	
	// Specific to this model
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVal3 &) override;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;
	using MetricElementType = NormType::HookeTensorType;
	
	// Generic
	typedef typename Traits::template DataSource<MetricElementType> MetricType;
	std::unique_ptr<MetricType> pMetric;
	ParamDefault param;
	
	ScalarType cosAngleMin = 0.5; // Refinement criterion
	virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(const PointType & p) const override;
	virtual DistanceGuess GetGuess(const IndexType & index) const override {return GetNorm(index);}
private:
	std::vector<OffsetType> tmp_stencil;
	NormType GetNorm(IndexCRef index) const; // Includes rescaling by h
};

#include "Implementation/Seismic2.hpp"


#endif /* Seismic2_h */
