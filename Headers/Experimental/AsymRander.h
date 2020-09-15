//
//  AsymRander.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 15/09/2020.
//

#ifndef AsymRander_h
#define AsymRander_h

#include "Base/Lagrangian2Stencil.h"
#include "Specializations/CommonTraits.h"
#include "JMM_CPPLibs/LinearAlgebra/AsymRanderNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"

#include "Seismic2.h" // Using the same generic stencil construction.

/// Traits and stencil for asym_rander_metrics, depending on a positive definite matrix and three vector fields.

struct TraitsAsymRander2 : TraitsBase<2> {
	typedef Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType> StencilType;
	typedef PeriodicGrid<TraitsAsymRander2> DomainType;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	using NormType = LinearAlgebra::AsymRanderNorm<ScalarType,2>;
	using DistanceGuess = NormType;
	using MetricElementType = NormType::AddableSourceType;
	using GridScalesType = PointType;
	
	static NormType MakeNorm(const MetricElementType & source, const GridScalesType & h){
		auto norm = NormType{source};
		norm.Rescale(h);
		return norm;
	}
	static constexpr GradientCachingStrategy useCache = GradientCachingStrategy::None;
};

using StencilAsymRander2 = StencilGenericLag2<TraitsAsymRander2>;


#endif /* AsymRander_h */
