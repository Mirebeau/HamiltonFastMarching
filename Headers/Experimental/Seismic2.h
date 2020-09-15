//
//  Seismic2.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 30/01/2019.
//

#ifndef Seismic2_h
#define Seismic2_h

#include <unordered_map>

#include "Base/Lagrangian2Stencil.h"
#include "Specializations/CommonTraits.h"

#include "JMM_CPPLibs/LinearAlgebra/SeismicNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/ComposedNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"

enum class GradientCachingStrategy {None,Local,Full};

struct TraitsSeismic2 : TraitsBase<2> {
	using StencilType = Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType>;
	using DomainType = PeriodicGrid<TraitsSeismic2>;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	using NormType = LinearAlgebra::SeismicNorm<ScalarType,Dimension>;
	using DistanceGuess = NormType;
	using MetricElementType = NormType::HookeTensorType;
	using GridScalesType = ScalarType;
	static NormType MakeNorm(const MetricElementType & m, GridScalesType h){
		return NormType{(1./square(h))*m};}
	
	static constexpr GradientCachingStrategy useCache = GradientCachingStrategy::Local;
};

struct TraitsSeismicTopographic2 : TraitsBase<2> {
	
	using StencilType = Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType>;
	using DomainType = PeriodicGrid<TraitsSeismicTopographic2>;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};

	using BaseNormType = LinearAlgebra::SeismicNorm<ScalarType,Dimension>;
	using TransformType = LinearAlgebra::Matrix<ScalarType,Dimension,Dimension>;
	using MetricElementType = LinearAlgebra::VectorPair<BaseNormType::HookeTensorType,TransformType>;
	using NormType = LinearAlgebra::ComposedNorm<BaseNormType,TransformType>;
	using DistanceGuess = NormType;
	
	using GridScalesType = PointType;
	static NormType MakeNorm(const MetricElementType & m, const GridScalesType h){
		TransformType a=m.second;
		for(int i=0; i<Dimension; ++i) {for(int j=0; j<Dimension; ++j) {a(i,j)*=h[j];}}
		return NormType{m.first, a};}
	
	static constexpr GradientCachingStrategy useCache = GradientCachingStrategy::Local;
};

template<typename TTraits>
struct StencilGenericLag2 final
: HamiltonFastMarching<TTraits>::StencilDataType {
	using Traits = TTraits;
	using HFM = HamiltonFastMarching<Traits>;
	using Superclass = typename HFM::StencilDataType;
	Redeclare13Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,
					IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,RecomputeType,
					DomainType,IndexDiff,PointType)
	Redeclare7Types(Traits,NormType,IndexType,StencilType,OffsetType,DistanceGuess,
					MetricElementType,GridScalesType)
	Redeclare1Type(Superclass,OffsetVal3)
	Redeclare1Constant(HFM,Dimension)
	
	// Specific to this model
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVal3 &) override;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;
	
	// Generic
	using MetricType = typename Traits::template DataSource<MetricElementType>;
	std::unique_ptr<MetricType> pMetric;
	using ParamType = std::conditional_t<std::is_same_v<GridScalesType,ScalarType>,
	ParamDefault,typename HFM::template _ParamDefault<2,void> >;
	ParamType param;
	
	ScalarType cosAngleMin = 0.5; // Refinement criterion
	virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(const PointType & p) const override;
	virtual DistanceGuess GetGuess(const IndexType & index) const override {return GetNorm(index);}
private:
	NormType GetNorm(IndexCRef index) const; // Includes rescaling by h
	
	std::vector<OffsetType> tmp_stencil;
	std::vector<VectorType> tmp_stencil_vec;
	std::vector<ScalarType> tmp_stencil_scal;
	template<size_t n> using Vec = typename NormType::template Vec<n>;
	
	// Tentative optimization : avoid recomputing gradients.
	// Actually counter productive : 50% more cpu time, cost of memory allocations exceeds gains.
	std::unordered_map<long,VectorType> vertexCache;
	std::unordered_multimap<DiscreteType,long> vertexCacheKeys;
	virtual void EraseCache(DiscreteType index) override final;
	static long hash(DiscreteType,OffsetType);
	template<typename Transform> struct MetricCaster;
	const GridScalesType & GridScales() const;
};

using StencilSeismic2 = StencilGenericLag2<TraitsSeismic2>;
using StencilSeismicTopographic2 = StencilGenericLag2<TraitsSeismicTopographic2>;

#include "Implementation/Seismic2.hpp"


#endif /* Seismic2_h */
