//
//  QuadLinLag2.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef QuadLinLag2_h
#define QuadLinLag2_h

#include "JMM_CPPLibs/Macros/DependentFalse.h"

#include "Base/Lagrangian2Stencil.h"
#include "JMM_CPPLibs/LinearAlgebra/RanderNorm.h"
#include "JMM_CPPLibs/LinearAlgebra/AsymmetricQuadraticNorm.h"

#include "Specializations/CommonTraits.h"
#include "JMM_CPPLibs/LinearAlgebra/HopfLaxMinimize.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"

/// Traits and stencil for metrics built of a quadratic and a linear part

struct TraitsRanderLag2 : TraitsBase<2> {
	typedef Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType> StencilType;
	typedef PeriodicGrid<TraitsRanderLag2> DomainType;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	typedef LinearAlgebra::RanderNorm<ScalarType,2> NormType;
	typedef LinearAlgebra::RanderNorm<ScalarType,1> NormType1;
	typedef NormType DistanceGuess;
};

struct TraitsAsymmetricQuadraticLag2 : TraitsBase<2> {
	typedef Lagrangian2Stencil<OffsetType,ScalarType,DiscreteType> StencilType;
	typedef PeriodicGrid<TraitsAsymmetricQuadraticLag2> DomainType;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};
	
	typedef LinearAlgebra::AsymmetricQuadraticNorm<ScalarType,2> NormType;
	typedef LinearAlgebra::AsymmetricQuadraticNorm<ScalarType,1> NormType1;
	typedef NormType DistanceGuess;
};

template<typename T>
struct StencilQuadLinLag2 final
: HamiltonFastMarching<T>::StencilDataType {
	typedef HamiltonFastMarching<T> HFM;
	typedef typename HFM::StencilDataType Superclass;
	Redeclare14Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,RecomputeType,Traits,PointType,
					 IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,DomainType,IndexDiff)
	Redeclare6Types(Traits,NormType,NormType1,IndexType,StencilType,OffsetType,DistanceGuess)
	Redeclare1Type(Superclass,OffsetVal3)
	Redeclare1Constant(HFM,Dimension)
	
	// Specific to this model
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef, const OffsetVal3 &) override;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) override;
	
	typedef typename NormType::SymmetricMatrixType SymmetricMatrixType;
	typedef LinearAlgebra::VectorPair<SymmetricMatrixType,VectorType> MetricElementType; // Allows averaging
	
	// Generic
	typedef typename Traits::template DataSource<MetricElementType> MetricType;
	std::unique_ptr<MetricType> pMetric;
	using ParamType = typename HFM::template _ParamDefault<2,void>; // Distinct scale on each axis
	ParamType param;

	virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(const PointType &) const override;
	virtual DistanceGuess GetGuess(const IndexType & index) const override {return GetNorm(index);}

	ScalarType cosAngleMin=0.5;
private:
	NormType GetNorm(IndexCRef) const;
	// Refinement near walls
	ScalarType wallBoundaryAngularResolution = 0.2;
	typename Traits::template Array<bool,Dimension> walls;
	const DomainType * pDom = nullptr;
	bool OnWallBoundary(IndexCRef) const;
	
	SymmetricMatrixType gridScales;
	NormType Rescale(const MetricElementType &) const;
	
	// Used by optimized Stern-Brocot algorithm.
	std::vector<OffsetType> tmp_stencil;
	std::vector<VectorType> tmp_stencil_vec;
	std::vector<ScalarType> tmp_stencil_scal;
	template<typename Norm, bool b> struct MetricCaster;
	template<bool b> void SetMetricCaster(HFMI *);
};

using StencilRander2 = StencilQuadLinLag2<TraitsRanderLag2>;
using StencilAsymmetricQuadratic2 = StencilQuadLinLag2<TraitsAsymmetricQuadraticLag2>;

#include "Implementation/QuadLinLag2.hxx"
#endif /* QuadLinLag2_h */
