//
//  Seismic3.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef Seismic3_h
#define Seismic3_h


#include "Base/Lagrangian3Stencil.h"
#include "Specializations/CommonTraits.h"

#include "JMM_CPPLibs/LinearAlgebra/SeismicNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/ComposedNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"


struct TraitsSeismic3 : TraitsBase<3> {
	using StencilType = Lagrangian3Stencil<OffsetType,ScalarType,DiscreteType>;
	using DomainType = PeriodicGrid<TraitsSeismic3>;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	using NormType = LinearAlgebra::SeismicNorm<ScalarType,3>;
	using DistanceGuess = NormType;
	
	using MetricElementType = NormType::HookeTensorType;
	static NormType MakeNorm(const MetricElementType & m, ScalarType h){
		return NormType{(1./square(h))*m};}
};

struct TraitsSeismicTopographic3 : TraitsBase<3> {
	
	using StencilType = Lagrangian3Stencil<OffsetType,ScalarType,DiscreteType>;
	using DomainType = PeriodicGrid<TraitsSeismicTopographic3>;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	using BaseNormType = LinearAlgebra::SeismicNorm<ScalarType,3>;
	using TransformType = LinearAlgebra::TopographicTransform<VectorType>;
	using NormType = LinearAlgebra::ComposedNorm<BaseNormType,TransformType>;
	using DistanceGuess = NormType;
	
	using MetricElementType = LinearAlgebra::VectorPair<BaseNormType::HookeTensorType,VectorType>;
	static NormType MakeNorm(const MetricElementType & m, ScalarType h){
		return NormType{(1./square(h))*m.first, TransformType{m.second}};}
	
};


template<typename TTraits>
struct StencilGenericLag3 final
: HamiltonFastMarching<TTraits>::StencilDataType {
	using Traits = TTraits;
	typedef HamiltonFastMarching<Traits> HFM;
	typedef typename HFM::StencilDataType Superclass;
	
	Redeclare15Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,
					IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,RecomputeType,
					DomainType,IndexDiff,PointType,IndexType,OffsetType)
	Redeclare4Types(Traits,NormType,StencilType,DistanceGuess,MetricElementType)
	Redeclare1Type(Superclass,OffsetVals)
	Redeclare1Constant(HFM,Dimension)
	
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVals &) override;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;
	
	typedef typename Traits::template DataSource<MetricElementType> MetricType;
	std::unique_ptr<MetricType> pMetric;
	ParamDefault param;
	bool checkAcuteness = false; // This is TODO

	virtual void SetStencil(IndexCRef index, StencilType & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(const PointType &) const override;
	virtual DistanceGuess GetGuess(const IndexType & index) const override {return GetNorm(index);}
private:
	NormType GetNorm(IndexCRef index) const; // Includes rescaling by h
	template<size_t n> using Vec = typename NormType::template Vec<n>;

	
	// Tentative optimization : Caching data for faster computations
	// Gets about 30 performance gain on small test case, for a much increase memory usage.
	// Not sure if worth it.
	const bool useHopfLaxCache = false;
	std::map<long,VectorType> vertexCache;
	std::map<long,std::pair<VectorType,ScalarType> > edgeCache;
	virtual void EraseCache(DiscreteType index) override final;
	static long hash(DiscreteType,OffsetType);
	static std::pair<long,bool> hash(DiscreteType,OffsetType,OffsetType);
};

using StencilSeismic3 = StencilGenericLag3<TraitsSeismic3>;
using StencilSeismicTopographic3 = StencilGenericLag3<TraitsSeismicTopographic3>;

#include "Implementation/Seismic3.hpp"


#endif /* Seismic3_h */
