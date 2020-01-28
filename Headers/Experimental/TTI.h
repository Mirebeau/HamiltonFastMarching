//
//  TTI.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 02/01/2020.
//

#ifndef TTI_h
#define TTI_h

#include "Base/Lagrangian2Stencil.h"
#include "Specializations/CommonTraits.h"

#include "JMM_CPPLibs/LinearAlgebra/TTINorm.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorPairType.h"
#include "JMM_CPPLibs/LinearAlgebra/BasisReduction.h"
#include "JMM_CPPLibs/LinearAlgebra/DifferentiationType.h"
#include "JMM_CPPLibs/LinearAlgebra/AD2.h"

template <int VDimension>
struct TraitsTTI : TraitsBase<VDimension> {
	using Superclass = TraitsBase<VDimension>;
	Redeclare3Types(Superclass,OffsetType,ScalarType,DiscreteType)
	struct StencilType;
	using DomainType = PeriodicGrid<TraitsTTI<VDimension> >;
	struct DifferenceType {static const int multSize = -1; struct MultiplierType {};};

	using NormType = LinearAlgebra::TTINorm<ScalarType,VDimension>;
	using DistanceGuess = NormType;
	Redeclare3Types(NormType,LinearType,QuadraticType,TransformType)
	using AlgebraicType = LinearAlgebra::VectorPair<LinearType,QuadraticType>;
	using MetricElementType = LinearAlgebra::VectorPair<AlgebraicType,TransformType>;
};

template<int VD> struct TraitsTTI<VD>::StencilType
: Lagrangian2Stencil<OffsetType, ScalarType, DiscreteType> {
	DiscreteType NSectors() const {assert(false);return 0;}
	static constexpr int Dimension = VD;
	static constexpr int nActiveNeigh = (Dimension*(Dimension+1))/2;
	using CommonStencilType = CommonStencil<OffsetType,ScalarType,nActiveNeigh>;
	Redeclare3Types(CommonStencilType,DiscreteFlowElement,DiscreteFlowType,RecomputeType);
	
    struct ActiveNeighFlagType {
		using SectorIndexType = void;
		ActiveNeighFlagType():t_(std::numeric_limits<ScalarType>::infinity()){};
        bool none() const {return t_==std::numeric_limits<ScalarType>::infinity();}

		ActiveNeighFlagType(ScalarType tActive, bool interior)
		:t_(interior ? tActive : -tActive){};
		bool isInterior(){return t_<0;}
		ScalarType tActive(){return std::abs(t_);}
		
		explicit ActiveNeighFlagType(ScalarType a):t_(a){}
        explicit operator ScalarType() const {return t_;}
	protected:
		ScalarType t_;
    };
};

template<int VD, typename Dummy> struct ExtraTraits<TraitsTTI<VD>,Dummy> {
	static constexpr bool hasBundle = false;
	static constexpr StencilStoragePolicy policy = StencilStoragePolicy::Lag2;
	static constexpr bool strictlyCausal = false;
	static constexpr bool factorFirstPass = true;
};

template<int VDimension>
struct StencilTTI final
: HamiltonFastMarching<TraitsTTI<VDimension> >::StencilDataType {
	using Traits = TraitsTTI<VDimension>;
	using HFM = HamiltonFastMarching<Traits>;
	using Superclass = typename HFM::StencilDataType;
	Redeclare14Types(HFM,IndexCRef,FullIndexCRef,OffsetType,OffsetCRef,ScalarType,
					 ActiveNeighFlagType,DiscreteType,DistanceGuess,PointType,IndexType,
					 ParamInterface,HFMI,RecomputeType,
					 DiscreteFlowType)
	Redeclare2Types(Traits,NormType,MetricElementType)
	Redeclare6Types(NormType,TransformType,SymmetricMatrixType,DiscreteVectorType,
					SellingPath,ReductionType,NeighborValuesType);
	Redeclare1Constant(HFM,Dimension)
	
	// Member fields
	using ParamType = typename HFM::template _ParamDefault<2,void>;
	ParamType param;
	virtual const ParamInterface & Param() const override {return param;}
	using MetricType = typename Traits::template DataSource<MetricElementType>;
	std::unique_ptr<MetricType> pMetric;
	virtual void Setup(HFMI *) override;

	// Specific to this model
	virtual void SetNeighbors(IndexCRef, std::vector<OffsetType> &) override;

	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType,
									 ActiveNeighFlagType &) override final;
	virtual RecomputeType
	_HopfLaxRecompute(IndexCRef, ActiveNeighFlagType, DiscreteFlowType &) override final;

	virtual DistanceGuess GetGuess(const PointType &) const override;
	virtual DistanceGuess GetGuess(const IndexType &) const override;
protected:
	NormType MakeNorm(const MetricElementType &) const;
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,
	const typename Superclass::OffsetVal3 &) override final {assert(false); return {};}
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &)
	override final {assert(false); return RecomputeType{};}


};

#include "Implementation/TTI.hpp"

#endif /* TTI_h */
