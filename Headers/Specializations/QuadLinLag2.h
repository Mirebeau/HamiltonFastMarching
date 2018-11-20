//
//  QuadLinLag2.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef QuadLinLag2_h
#define QuadLinLag2_h

#include "Base/Lagrangian2Stencil.h"

/// Stencil for metrics built of a quadratic and a linear part
template<typename T>
struct StencilQuadLinLag2 : HamiltonFastMarching<T>::StencilDataType {
	typedef HamiltonFastMarching<T> HFM;
	typedef typename HFM::StencilDataType Superclass;
	Redeclare6Types(HFM,ParamDefault,ParamInterface,HFMI,DiscreteFlowType,RecomputeType,Traits)
	Redeclare7Types(HFM,IndexCRef,VectorType,ScalarType,DiscreteType,OffsetCRef,DomainType,IndexDiff)
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
	ParamDefault param;
	
	virtual void SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) override;
	virtual const ParamInterface & Param() const override {return param;}
	virtual void Setup(HFMI *) override;
	virtual DistanceGuess GetGuess(IndexCRef index) const override;
private:
	std::forward_list<OffsetType> l;
	NormType GetNorm(IndexCRef index) const {
		const MetricElementType data = (*pMetric)(index);
		return NormType{data.first,data.second};}
	
	// Refinement near walls
	ScalarType wallBoundaryAngularResolution = 0.2;
	typename Traits::template Array<bool,Dimension> walls;
	/*    typedef typename HFM::template DataSource<bool> BoolField;
	 std::unique_ptr<BoolField> pWalls = nullptr;*/
	const DomainType * pDom = nullptr;
	bool OnWallBoundary(IndexCRef) const;
};


typedef StencilQuadLinLag2<TraitsRanderLag2> StencilRanderLag2;
typedef StencilQuadLinLag2<TraitsAsymmetricQuadraticLag2> StencilAsymmetricQuadraticLag2;

#include "QuadLinLag2.hxx"
#endif /* QuadLinLag2_h */
