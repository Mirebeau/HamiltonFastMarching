//
//  StaticFactoring.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 17/09/2020.
//

#ifndef StaticFactoring_h
#define StaticFactoring_h

/**
 This file implements static additive factoring for eikonal equations, a technique used to mitigate the impact of the solution
 singularities on the accuracy of the numerical scheme. Static factorization is only applied at the source points only, and the
 other singularities (obstacle corners) are not dealt with.
 
 Static factorization is intended as a replacement for dynamic factorization, a more complicated technique intended for
 dealing (usually poorly) with obstacles, and which was previously implemented in an over-enthusiastic moment.
  */

template<typename T>
struct StaticFactoring {
	typedef T Traits;
    typedef HamiltonFastMarching<Traits> HFM;
    Redeclare16Types(HFM,IndexCRef,FullIndexCRef,DiscreteFlowType,OffsetType,OffsetCRef,
					 ScalarType,DiscreteType,IndexType,IndexDiff,DomainTransformType,
					 VectorType,HFMI,PointType,DistanceGuess,DiscreteFlowElement,DomainType)
	Redeclare1Constant(HFM,Dimension)
    template<typename E, size_t n> using Array = typename Traits::template Array<E,n>;
	
	Array<bool,Dimension> mask;
	Array<ScalarType,Dimension> values;
	Array<VectorType,Dimension> gradients;

	/// Wether the arrays cover only a subdomain, and the indexShift in that case.
	bool subdomain; IndexDiff indexShift;
	bool Disabled() const {return values.empty();}
	bool recompute_default;

	/// Returns wether recomputation is needed in the FMM, before freezing the point
	bool NeedsRecompute(IndexType ind) const;
	
	/// Returns wether factorization is active at given point, and sets the index if appropriate
	bool SetIndex(IndexCRef ind) {active=_SetIndex(ind); return active;}
	bool _SetIndex(IndexCRef ind);
	
	/** Correction to apply to a finite difference at the given position, with the given offset, and order.
	Expects :
	   - the factorization applies to the current index
	   - ind + i * offset is in the domain, for all 0<=i <=order.
	 */
	ScalarType Correction(const OffsetType & off, int order) const;
	
    void Setup(HFMI *);
	ScalarType seedRadius;
protected:
	IndexType currentIndex; // Shifted if needed
	bool active=false; // Wether factorization applies at the current index
	const HFM * pFM; // Reference domain, for periodization
	void ComputeFactor(HFMI*); // Sets the factor from the metric
	void ImportFactor(HFMI*);  // Imports a given factor
	void SetSeeds(HFMI*);
	PointType seed;
};

#include "Implementation/StaticFactoring.hxx"

#endif /* StaticFactoring_h */
