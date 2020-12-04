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


template<typename T> void StaticFactoring<T>::
Setup(HFMI * that){
	IO & io = that->io;
	pFM = that->pFM.get();
	const HFM & fm = *pFM;
	recompute_default = fm.order>1 || ! HFM::factorFirstPass;

	const bool
	import = io.HasField("factoringValues"),
	compute = io.HasField("factoringRadius");
	if(import && compute){
		ExceptionMacro("StaticFactoring error : factoringValues (import)"
		" and factoringRadius (compute) both specified.")}
	if(import){ImportFactor(that);}
	SetSeeds(that);
	if(compute){ComputeFactor(that);}
}

template<typename T> void StaticFactoring<T>::
ImportFactor(HFMI * that){
	IO & io = that->io;
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
	const IndexType dims= subdomain ? values.dims : pFM->values.dims;
	if(values.dims!=dims || gradients.dims!=dims || mask.dims!=dims){
		ExceptionMacro("Static factoring error : inconsistent data dimension.")}
		
	// Rescale the gradients. The gradient is a co-vector, which must be adimensionized.
	// This amounts to re-dimensioning a vector.
	const auto & param = pFM->stencilData.Param();
	for(VectorType & g : gradients) {g = param.ReDim(g);}
}

template<typename T> void StaticFactoring<T>::
ComputeFactor(HFMI * that){
	// This function sets the factor if it is not provided by the user
	const HFM & fm = *pFM;
	const auto & dom = fm.dom;
	const auto & stencil = fm.stencilData;
	auto & io = that->io;
	
	
	io.SetHelp("factoringRadius", "Radius, in pixels, where to use source factorization");
	const DiscreteType factoringRadius = io.template Get<ScalarType>("factoringRadius",10);
	
	IndexType dims;
	if(factoringRadius<0){
		subdomain=false;
		dims = fm.values.dims;
		indexShift=IndexDiff::Constant(0);
	} else {
		subdomain=true;
		IndexType bottom,top;
		const IndexType dimax = fm.values.dims;
		for(int i=0; i<Dimension; ++i){
			bottom[i] = std::max(0,(DiscreteType)floor(seed[i]-factoringRadius));
			top[i]    = std::min(dimax[i]-1,(DiscreteType)ceil(seed[i]+factoringRadius));
		}
		indexShift = IndexDiff::FromOrigin(bottom);
		dims = top-bottom;
	}
	
	// Reshape the arrays
	const ScalarType NaN = std::numeric_limits<ScalarType>::quiet_NaN();
	const DiscreteType size = dims.Product();
	values.dims = dims;    values.resize(size,NaN);
	gradients.dims = dims; gradients.resize(size,VectorType::Constant(0.));
	mask.dims = dims;      mask.resize(size,true);

	io.SetHelp("factoringPointChoice", "Choice of source factorization (Key,Both)");
	const std::string choice
	= io.GetString("factoringPointChoice", fm.order<3 ? "Seed" : "Both");

	const auto & seedDist = stencil.GetGuess(seed);
		
	// Note : asymmetric norms in Python and c++ are defined with opposite asymmetry,
	// aka norm_python(v) = norm_c++(-v). This is taken into account below.
	
	if(choice=="Seed" || choice=="Key"){ // "Seed" and "Key" are synonyms here
		// Factorization dist_{x_0}(x-x_0), where x_0 is the seed
		//TODO : it would be preferable to do this loop in parallel. However, a convenient
		//implementation requires C++ 20, which is not widely supported yet.
		// https://stackoverflow.com/a/52834495/12508258
		for(int i=0; i<size; ++i){
			const PointType p = dom.PointFromIndex(values.Convert(i) + indexShift);
			const VectorType v = p-seed;
			// " - grad( - v) " is because of the different orientation convention between c++ and Python, for asymmetric norms
			const VectorType g = v.IsNull()? VectorType::Constant(0.):
				-seedDist.Gradient(-v);
			gradients[i] = g;
			values[i] = g.ScalarProduct(v);
		}
	} else if(choice=="Both"){
		// Factorization (dist_{x_0}(x-x_0) + dist_{x}(x-x_0)) / 2
		// Not parallelisable as is
		
		// Issue : the bilinear interpolation of the metric at the seed defeats the
		// purpose of this factorization choice, which is third order accuracy.
		// Also the use of centered finite differences for the distance gradient w.r.t
		// position may not be accurate enough for third order accuracy.
		// Therefore, using 'Both' choice in the c++ code is discouraged at this point,
		// and this call should be intercepted by the Python interface, which uses higher
		// interpolation methods (third order splines). See dictIn.SetFactor
		
		WarnMsg() << "factoringPointChoice = 'Both' yields a factor that is "
		"best computed in Python, via dictIn.SetFactor, "
		"which uses higher order interpolation methods.";
		
		for(int pi=0; pi<size; ++pi){
			const IndexType pIndex = values.Convert(pi) ;
			const PointType p = dom.PointFromIndex(pIndex + indexShift);
			const VectorType v = p-seed;
			const auto & pDist = stencil.GetGuess(pIndex);
			const VectorType seedG = v.IsNull()?VectorType::Constant(0.):
				-seedDist.Gradient(-v);
			const VectorType pG = v.IsNull()?VectorType::Constant(0.):
				-pDist.Gradient(-v);
			const VectorType g = 0.5*(pG+seedG);
			values[pi] = g.ScalarProduct(v);
			
			// We need to differentiate
			// p,v -> (dist_{x_0}(v) + dist_{p}(v)) / 2
			// w.r.t p,v, defined as p=x, v = x-x_0
			
			// Differentiating w.r.t v
			gradients[pi]+=g;
			
			// Differentiating w.r.t p using
			// centered finite differences in interior,
			// and first order finite diff along boundary
			for(int k=0; k<Dimension; ++k){
				for(int eps=-1; eps<=1; eps+=2){
					IndexType qIndex = pIndex;
					qIndex[k]+=eps;
					if(mask.InRange(qIndex)){
						const PointType q = dom.PointFromIndex(qIndex + indexShift);
						ScalarType qVal = pDist.Norm(-(q-seed))/2.;
						const DiscreteType qi = gradients.Convert(qIndex);
						qIndex[k]+=eps;
						if(mask.InRange(qIndex)) qVal/=2; // True iff centered finite diff
						gradients[qi][k] -= eps * qVal;
					} else {
						gradients[pi][k] += eps * pG.ScalarProduct(v)/2.; // Non centered diff
					}
				} // for eps
			} // for k
		} // for pi
	} else {
		ExceptionMacro("StaticFactoring error : Unsupported factoringPointChoice "
					   << choice << ".");
	}
	if(io.template Get<ScalarType>("exportFactoring",0.)){
		io.template SetArray<ScalarType,Dimension>("factoringValues",values);
		io.template SetArray<VectorType,Dimension>("factoringGradients",gradients);
		if(subdomain) io.template Set<VectorType>("factoringIndexShift",
												  VectorType::CastCoordinates(indexShift));
	}
}


template<typename T> void StaticFactoring<T>::
SetSeeds(HFMI * that) {
	auto & io = that->io;
	const auto & stencil = pFM->stencilData;
	std::vector<PointType> seedPoints;
	std::vector<ScalarType> seedValues;
	
	const bool compute = io.HasField("factoringRadius");
	const bool disabled = Disabled() && ! compute; // Only ImportFactor run
	
	io.SetHelp("seeds",
			   "Type: vector of points. (Input)\n"
			   "Usage: sets the points from which the front propagation starts.");
	
	if(io.HasField("seeds")){
		seedPoints = io.template GetVector<PointType>("seeds");
		for(auto & p : seedPoints) p=stencil.Param().ADim(p);
		if(compute){
			if(seedPoints.size()==1){seed=seedPoints[0];}
			else {ExceptionMacro("Exactly one seed point must be specified for computing "
								 "the factor, found "<< seedPoints.size() << ".");}
		}
		
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
		
		seedRadius = io.template Get<ScalarType>("seedRadius",disabled ? 0. : 2.);
		const bool stencilStep = seedRadius<0;
		seedRadius = std::abs(seedRadius);
		
		if(seedRadius>0){
			std::vector<PointType> newPoints;
			std::vector<ScalarType> newValues;
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
		if(compute && ! uPoints.empty()){
			ExceptionMacro("Factor computation incompatible with unoriented seed points");}
		
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
			SpecializationsDefault::PadAdimEquiv(that,uPoints[i],equiv);
			seedPoints.insert(seedPoints.end(),equiv.begin(),equiv.end());
			seedValues.resize(seedValues.size()+equiv.size(),uValues[i]);
		}
	}
	
	for(int i=0; i<seedPoints.size(); ++i){
		auto seedIndex=pFM->dom.IndexFromPoint(seedPoints[i]);
		if(!pFM->dom.PeriodizeNoBase(seedIndex).IsValid()){
			WarnMsg() << "Error : seed " << stencil.Param().ReDim(seedPoints[i]) << " is out of range.\n";
			continue;}
		that->pFM->seeds.insert({seedIndex,seedValues[i]});
	}
	if(pFM->seeds.empty() && !seedPoints.empty())
		ExceptionMacro("Error : seeds incorrectly set");
}

#endif /* StaticFactoring_hxx */
