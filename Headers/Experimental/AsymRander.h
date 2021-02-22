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

/*
// This comment shows how one may enfore an upper bound on the radius of the semi-Lagrangian
 stencils during refinement (at the cost of some non-causality in strongly anisotropic cases)
 
struct StencilAsymRander2 : StencilGenericLag2<TraitsAsymRander2> {
	virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
	const int squaredStencilRadiusBound=10; // Bound on the squared stencil radius
};

void StencilAsymRander2::SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil){
	const NormType & norm = GetNorm(index);
	assert(tmp_stencil.empty());
	tmp_stencil.insert(tmp_stencil.end(),{OffsetType(1,0),OffsetType(0,-1),OffsetType(-1,0),OffsetType(0,1)});
	
	auto pred = [&norm,this](OffsetCRef u, OffsetCRef v) -> bool {
		// Do not refine if that would exceed the enforced bound on stencil radius
		if ((u+v).SquaredNorm()>squaredStencilRadiusBound) return false;
		// Refine if the cos(angle) is above desired bound (in particular if obtuse)
		const auto u_=VectorType::CastCoordinates(u), v_=VectorType::CastCoordinates(v);
		return CosAngle(norm,u_,v_) >= this->cosAngleMin;};
	
	SternBrocotRefine(pred, stencil, tmp_stencil);
}
*/

#endif /* AsymRander_h */
