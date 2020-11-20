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
    Redeclare15Types(HFM,IndexCRef,FullIndexCRef,DiscreteFlowType,OffsetType,
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

	/// Returns wether recomputation is needed in the FMM
	bool NeedsRecompute(IndexType ind) const {
		if(Disabled()) return false;
		if(subdomain){
			ind-=indexShift;
			if( ! mask.InRange(ind) ) {return false;}
		}
		return mask(ind) && recompute_default;
	}
	
	/// Returns wether factorization is active at given point, and sets the index if appropriate
	bool SetIndex(IndexCRef ind) {active=_SetIndex(ind); return active;}
	bool _SetIndex(IndexCRef ind) {
		if(Disabled()) return false; // No factorization applied
		currentIndex = ind;
		if(subdomain){
			currentIndex -= indexShift;
			if( ! mask.InRange(currentIndex) ) {return false;}
		}
		return mask(currentIndex);
	}
	
	/// Correction to apply to a finite difference at the given position, with the given offset, and order.
	ScalarType Correction(const OffsetType & off, int order) const {
		if(!active) return 0.;
	/* Input guarantees :
		- the factorization applies to the current index
		- ind + i * offset is in the domain, for all 0<=i <=order.
	 */
		const DomainType & dom = pFM->dom;
		const IndexType index  = currentIndex;
		const IndexDiff offset = IndexDiff::CastCoordinates(off);
		IndexType ind_ = index + order*offset; // Last index of Finite Difference
		if(subdomain && !(mask.InRange(ind_))){return 0;} // Out of domain FD
		else {dom.Periodize(ind_,index);}
		
		// Note : the gradient and offset associated to same position ind, so there is no
		// need for transform.PullVector etc
		const ScalarType deriv =
		gradients(index).ScalarProduct(VectorType::CastCoordinates(offset));
		
		if(order==1) {return deriv + (values(index) - values(ind_));}
		
		IndexType ind1 = index + offset;
		if(!subdomain) dom.Periodize(ind1,index);
		if(order==2) {return deriv
			+ (1.5*values(index) - 2*values(ind1) + 0.5*values(ind_));}
		
		IndexType ind2 = index + 2*offset;
		if(!subdomain) dom.Periodize(ind2,index);
		assert(order==3); // Order must be 1, 2 or 3.
		return deriv
		+ ((11./6.)*values(index)-3.*values(ind1)+1.5*values(ind2)-(1./3.)*values(ind_));
	}
	
    bool Setup(HFMI *);
protected:
	IndexType currentIndex; // Shifted if needed
	bool active=false; // Wether factorization applies at the current index
	const HFM * pFM; // Reference domain, for periodization
//	void SetFactor(); // Sets the factor from the metric
};


template<typename T> bool
StaticFactoring<T>::Setup(HFMI * that){
	IO & io = that->io;
	pFM = that->pFM.get();
	const HFM & fm = *pFM;
	recompute_default = fm.order>1 || ! HFM::factorFirstPass;

	// Import data (mask, values, gradient, indexShift), which is computed off site.
	
	if( ! io.HasField("factoringValues") ){ return false;}
	io.SetHelp("factoringValues","Values used in additive source factorization");
	values = io.GetArray<ScalarType,Dimension> ("factoringValues");
	
	io.SetHelp("factoringGradients","Gradients used in additive source factorization");
	if( ! io.HasField("factoringGradients") )
		ExceptionMacro("StaticFactoring error : factoringGradients are not provided");
	gradients = io.GetArray<VectorType,Dimension>("factoringGradients");

	io.SetHelp("factoringMask","Mask used in additive source factorization");
	if( ! io.HasField("factoringMask") ){
		mask.dims = values.dims; mask.resize(values.size(),true);
	} else {
		mask = io.GetArray<ScalarType,Dimension>("factoringMask").template Cast<bool>();
	}

	io.SetHelp("factoringIndexShift","Shift used for the source factorization data");
	subdomain = io.HasField("factoringIndexShift");
	indexShift = !subdomain ? IndexDiff::Constant(0) :
		IndexDiff::CastCoordinates(io.Get<VectorType>("factoringIndexShift"));
	
	// Check data
	const IndexType dims= subdomain ? values.dims : fm.values.dims;
	if(values.dims!=dims || gradients.dims!=dims || mask.dims!=dims){
		ExceptionMacro("Static factoring error : inconsistent data dimension.")}
		
	// Rescale the gradients. The gradient is a co-vector, which must be adimensionized.
	// This amounts to re-dimensioning a vector.
	const auto & param = fm.stencilData.Param();
	for(VectorType & g : gradients) {g = param.ReDim(g);}
	return true;
}

/*
template<typename T> void
StaticFactoring<T>::SetFactor(){
	// This function sets the factor if it is not provided by the user
	const HFM & fm = *pFM;
	const auto & dom = fm.dom;
	const auto & stencil = fm.stencil;
	
	// Reshape the arrays
	const ScalarType NaN = std::numeric_limits<ScalarType>::quiet_Nan();
	const IndexType dims=fm.values.dims;
	const DiscreteType size = fm.values.size();
	values.dims = dims;    values.resize(size,NaN);
	gradients.dims = dims; gradients.resize(size,VectorType::Constant(NaN));
	mask.dims = dims;      mask.resize(size,true);

	// Basic factorization : for each point, evaluate the gradient, deduce the norm
	if(fm.seeds.size()!=1)
	const PointType seed; assert(false); // TODO
	const auto & distSeed = stencil.GetGuess(seed);
	
	//TODO : it would be preferable to do this in parallel. However, a convenient
	//implementation requires C++ 20, which is not widely supported yet.
	// https://stackoverflow.com/a/52834495/12508258 
	for(int i=0; i<size; ++i){
		const PointType p = dom.PointFromIndex(values.Convert(i));
		const VectorType v = p-seed;
		const VectorType g = distSeed.Gradient(v);
		gradients[i] = g;
		values[i] = g.ScalarProduct(v);
	}
	
	assert(false);
	for(int i=0; i<size; ++i){
		const PointType p = dom.PointFromIndex(values.Convert(i));
		const VectorType v = p-seed;
		const auto distp =
	}
	

}
*/
#endif /* StaticFactoring_h */
