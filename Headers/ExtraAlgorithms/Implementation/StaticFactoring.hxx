//
//  StaticFactoring.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2020.
//

#ifndef StaticFactoring_hxx
#define StaticFactoring_hxx

template<typename T> bool StaticFactoring<T>::
NeedsRecompute(IndexType ind) const {
	if(Disabled()) return false;
	if(subdomain){
		ind-=indexShift;
		if( ! mask.InRange(ind) ) {return false;}
	}
	return mask(ind) && recompute_default;
}

template<typename T> bool StaticFactoring<T>::
_SetIndex(IndexCRef ind) {
	if(Disabled()) return false; // No factorization applied
	currentIndex = ind;
	if(subdomain){
		currentIndex -= indexShift;
		if( ! mask.InRange(currentIndex) ) {return false;}
	}
	return mask(currentIndex);
}


template<typename T> auto StaticFactoring<T>::
Correction(const OffsetType & off, int order) const -> ScalarType {
	if(!active) return 0.;
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


template<typename T> bool StaticFactoring<T>::
Setup(HFMI * that){
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
template<typename T> void StaticFactoring<T>::
SetFactor(){
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

template<typename T> void StaticFactoring<T>::
SetSeeds(HFMI * that) {
	const auto & io = that->io;
	const auto & stencil = pFM->stencil;
	std::vector<PointType> seedPoints;
	std::vector<ScalarType> seedValues;
	
	io.SetHelp("seeds",
			   "Type: vector of points. (Input)\n"
			   "Usage: sets the points from which the front propagation starts.");
	if(io.HasField("seeds")){
		seedPoints = io.template GetVector<PointType>("seeds");
		for(auto & p : seedPoints) p=stencil.Param().ADim(p);
		
		io.SetHelp("seedValues",
				   "Type: vector of scalars, same length as seeds. (Input)\n"
				   "Usage: sets the values at which the front propagation starts.");
		if(io.HasField("seedValues")) {
			seedValues = io.template GetVector<ScalarType>("seedValues");
			if(seedValues.size()!=seedPoints.size())
				ExceptionMacro("Error : Inconsistent size of seedValues.")
				} else {
					seedValues.resize(seedPoints.size(),0.);
				}
		
		io.SetHelp("seedRadius",
				   R"(Type: scalar. Unit: pixels. (Input)
				   Usage: The seedRadius option allows to spread seeds on the grid, instead of
				   rounding the given position to the nearest grid point (default).
				   The grid points within the given radius in pixels of the seed position are
				   set as seeds for the grid propagation, with appropriate values.
				   Special case: If radius is negative, then we also include those points
				   accessible in one step from the reverse stencil.)");
		
		// Setting default issue : factoring is not yet set.
		// Thus cannot test pFM->factoring.method == FactoringMethod::None
		bool factors = io.HasField("factoringValues");
		if(io.HasField("factoringMethod")) {
			factors = io.GetString("factoringMethod")!="None";}
		seedRadius = io.template Get<ScalarType>("seedRadius",factors ? 2. : seedRadius);
		const bool stencilStep = seedRadius<0;
		seedRadius = std::abs(seedRadius);
		
		if(seedRadius>0){
			std::vector<PointType> newPoints;
			std::vector<ScalarType> newValues;
			//				std::vector<VectorType> newGradients;
			// Avoid repetition of spreaded points, for a given seed point
			std::set<IndexType> indices;
			
			// Enumerate candidate offsets in a box
			Array<int, Dimension> arr;
			const DiscreteType seedRadiusBd = 1+(DiscreteType)std::ceil(seedRadius);
			arr.dims.fill(1+2*seedRadiusBd);
			IndexType arrCenter; arrCenter.fill(seedRadiusBd);
			auto arrOffset = [&](DiscreteType i){return arr.Convert(i)-arrCenter;};
			
			const auto & dom = pFM->dom;
			for(size_t i=0; i<seedPoints.size(); ++i){
				indices.clear();
				const PointType & p = seedPoints[i];
				const ScalarType value = seedValues[i];
				const auto & distp = stencil.GetGuess(p);
				constexpr ScalarType inf=std::numeric_limits<ScalarType>::infinity();
				
				auto insert = [&](IndexType index,ScalarType radius){
					const PointType q = dom.PointFromIndex(index);
					const VectorType v=p-q; // toward seed
					if(v.Norm()>=radius) return;
					if(!indices.insert(index).second) return; // Already seen
					IndexType perIndex = index;
					if(!pFM->dom.PeriodizeNoBase(perIndex).IsValid()) return; // out
					const auto & distq = stencil.GetGuess(perIndex);
					
					newPoints.push_back(q);
					newValues.push_back(value + 0.5*(distp.Norm(v)+distq.Norm(v)));
				}; // insert
				
				const IndexType pIndex = dom.IndexFromPoint(p);
				insert(pIndex,inf);
				
				for(DiscreteType i=0; i<arr.dims.Product(); ++i){
					insert(pIndex+arrOffset(i),seedRadius);}
				
				if(!stencilStep) {continue;}
				std::set<IndexType> indices2 = indices; // Copy already inserted
				for(IndexCRef index : indices2){
					const auto rev =
					stencil.ReversedOffsets({index,pFM->values.Convert(index)});
					for(OffsetCRef offset : rev){
						insert(index+IndexDiff::CastCoordinates(offset),inf);
					}
				}
			}
			std::swap(seedPoints,newPoints);
			std::swap(seedValues,newValues);
			
			// Export generated data
			std::vector<PointType> seedPointsRedim;
			for(const PointType & p : seedPoints) {seedPointsRedim.push_back(stencil.Param().ReDim(p));}
			io.SetVector("spreadedSeeds", seedPointsRedim);
			io.SetVector("spreadedSeedValues", seedValues);
			//				io.SetVector("spreadedSeedGradients", newGradients);
		}
	}
	
	using SpecializationsDefault = typename HFMI::SpecializationsDefault;
	
	if(HFM::hasBundle){
		io.SetHelp("seeds_Unoriented",
				   R"(Type: vector of points in the *base horizontal domain*. (Input)
				   Usage: set the points from which the front propagation starts.
				   An unoriented seed is equivalent to a family of standard seeds,
				   with the same base position and each possible bundle direction.)");}
	if(HFM::hasBundle && io.HasField("seeds_Unoriented")){
		using UnorientedPointType = typename SpecializationsDefault::UnorientedPointType ;
		const auto uPoints = io.template GetVector<UnorientedPointType>("seeds_Unoriented");
		
		std::vector<ScalarType> uValues;
		if(io.HasField("seedValues_Unoriented")){
			uValues = io.template GetVector<ScalarType>("seedValues_Unoriented");
		} else {
			uValues.resize(uPoints.size(),0.);
		}
		io.SetHelp("seedValues_Unoriented",
				   R"(Type: vector of scalar values. (Input)
				   Usage: set the values at which the front propagation
				   tarts, for each member of seeds_Unoriented.)");
		if(uValues.size()!=uPoints.size()) ExceptionMacro("Error : Inconsistent size of seedValues_Unoriented.");
		
		std::vector<PointType> equiv;
		for(int i=0; i<uPoints.size(); ++i){
			SpecializationsDefault::PadAdimEquiv(this,uPoints[i],equiv);
			seedPoints.insert(seedPoints.end(),equiv.begin(),equiv.end());
			seedValues.resize(seedValues.size()+equiv.size(),uValues[i]);
		}
	}
	
	for(int i=0; i<seedPoints.size(); ++i){
		auto seedIndex=pFM->dom.IndexFromPoint(seedPoints[i]);
		if(!pFM->dom.PeriodizeNoBase(seedIndex).IsValid()){
			WarnMsg() << "Error : seed " << stencil.Param().ReDim(seedPoints[i]) << " is out of range.\n";
			continue;}
		pFM->seeds.insert({seedIndex,seedValues[i]});
	}
	if(pFM->seeds.empty() && !seedPoints.empty())
		ExceptionMacro("Error : seeds incorrectly set");
}

#endif /* StaticFactoring_hxx */
