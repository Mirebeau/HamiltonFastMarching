//
//  TTI.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 07/01/2020.
//

#ifndef TTI_hpp
#define TTI_hpp

template<int VD> void StencilTTI<VD>::
Setup(HFMI * that) {
	Superclass::Setup(that);
	param.Setup(that);
	pMetric = that->template GetField<MetricElementType>("metric",false);
	auto & io = that->io;
	io.SetHelp("nmix","Number of ellipsoids defining the envelope of the TTI metric. "
			   "Special value nmix=0 : newton-like optimization (default).");
	if(io.HasField("nmix")){nmix = io.template Get<ScalarType>("nmix");}
	if(nmix<0){ExceptionMacro("Expect non-negative nmix, found " << nmix);}
}

template<int VD> auto StencilTTI<VD>::
MakeNorm(const MetricElementType & m) const -> NormType {
	const auto & [alg,trans] = m;
	const auto & [linear,quadratic] = alg;
	NormType norm{linear,quadratic,trans};
	for(int i=0; i<Dimension; ++i){
		for(int j=0; j<Dimension; ++j){
			norm.transform(i,j)*=param.gridScales[j];
		}
	}
	return norm;
}

template<int VD> auto StencilTTI<VD>::
GetGuess(const PointType & p) const -> DistanceGuess {
	return MakeNorm(MapWeightedSum<MetricElementType>(*pMetric,
													  this->pFM->dom.Neighbors(p)));
}

template<int VD> auto StencilTTI<VD>::
GetGuess(const IndexType & index) const -> DistanceGuess {
	return MakeNorm( (*pMetric)(index) );
}

template<int VD> void StencilTTI<VD>::
SetNeighbors(IndexCRef index, std::vector<OffsetType> & neigh) {
	// Returns the stencil used by the Eulerian scheme.
	// Opposite offsets are added.
	// No check for redundancy.
	const NormType norm = GetGuess(index);
	const auto props = norm.Props();
	auto selling = norm.Selling(props.tMin);
	
	// Push the offsets.
	auto push_opp = [&neigh](const DiscreteVectorType & v){
		neigh.push_back(OffsetType::CastCoordinates(v));
		neigh.push_back(OffsetType::CastCoordinates(-v));
	};

	for(const auto & e : selling.offsets){push_opp(e);}
	
	// If riemannian, return.
	if(props.tMax==props.tMin) return;
	
	ScalarType t=props.tMin;
	while(true){
		const auto [i,j] = selling.NextStep(t);
		// If no Selling step before tMax, return.
		if(t>=props.tMax) break;
		selling.MakeStep({i,j});
		push_opp(selling.offsets[ReductionType::LinearizeIndices(i,j)]);
	}
}

template<int VD> auto StencilTTI<VD>::
HopfLaxUpdate_nmix(FullIndexCRef updated, OffsetCRef acceptedOffset, ScalarType value,
				   ActiveNeighFlagType & active) -> ScalarType {
	auto & fm = *(this->pFM);
	const ScalarType oldValue = fm.values[updated.linear];
	fm.template SetIndex<true,false>(updated.index); // useFactoring, smallCorrection

	const NormType norm = GetGuess(updated.index);
	const auto props = norm.Props();
	const bool mix_is_min = props.optimDirection==-1;
	const ScalarType inf = std::numeric_limits<ScalarType>::infinity();

	const TransformType & A = norm.transform.Inverse();
	using Sym = SymmetricMatrixType;
	const Sym D0 = Dimension==2 ? Sym::RankOneTensor(A.Row(0)) :
	Sym::RankOneTensor(A.Row(0)) + Sym::RankOneTensor(A.Row(1));
	const Sym D1 = Sym::RankOneTensor(A.Row(Dimension-1));
	
	ScalarType tOpt,
	solOpt = mix_is_min ? inf : -inf;
	
	for(int imix = 0; imix<nmix; ++imix){
		const ScalarType t =
		nmix==1 ? (props.tMin+props.tMax)/2 :
		props.tMin + (props.tMax-props.tMin)* imix/ScalarType(nmix-1);
		const Sym D = (1-t)*D0+t*D1;
		const ScalarType mult = norm.Multiplier(t);
		const auto & decomp = ReductionType::TensorDecomposition(D);

		// Get the neighbor values
		const int KKTDim = ReductionType::KKTDimension;
		std::array<ScalarType,KKTDim> val;
		for(int k=0; k<KKTDim; ++k){
			int orderp=1,orderm=1;
			const auto offset = OffsetType::CastCoordinates(decomp.offsets[k]);
			const ScalarType
			valp = fm.template GetNeighborValue<true,false,1>( offset,orderp),
			valm = fm.template GetNeighborValue<true,false,1>(-offset,orderm);
			val[k] = std::min(orderp ? valp : inf, orderm ? valm : inf);
		}
		
		// Sort and solve
		std::array<ScalarType,KKTDim> indices;
		for(int i=0; i<KKTDim; ++i) {indices[i]=i;}
		std::sort(indices.begin(),indices.end(),
				  [&val](int i,int j){return val[i]<val[j];});
		const ScalarType valMin=val[indices[0]];
		if(valMin==inf) continue;
			
		using std::sqrt; using std::max;
		ScalarType a(0),b(0),c(-mult),sol(inf);
		int r=0;
		for(; r<KKTDim; ++r){
			const int i = indices[r];
			// Optimization opportinity : first loop yiels v=0 (but w!=0)
			const ScalarType v = val[i] - valMin;
			if(v>=sol) break;
			
			const ScalarType w = decomp.weights[i], wv=w*v, wvv=wv*v;
			a+=w;
			b+=wv;
			c+=wvv;
			
			if(a<=mult*1e-10) {continue;}
			const ScalarType delta = b*b-a*c;
			assert(delta>=0.);
			sol = (b+sqrt(delta))/a;
		}
		sol+=valMin;
		
		// Extremize
		if(mix_is_min == (sol<solOpt)) {tOpt = t; solOpt=sol;}
	}
	
	const ScalarType newValue = solOpt;
	// newValue==-inf may arise if none of the neighbors are accepted, and mix is max
	if(newValue<oldValue && newValue!=-inf) {
		active = ActiveNeighFlagType(tOpt,true);
		return newValue;
	} else {
		return oldValue;
	}
}

template<int VD> auto StencilTTI<VD>::
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef acceptedOffset, ScalarType value,
			  ActiveNeighFlagType & active) -> ScalarType {
	if(nmix>0) {return HopfLaxUpdate_nmix(updated,acceptedOffset,value,active);}
	auto & fm = *(this->pFM);
	const ScalarType oldValue = fm.values[updated.linear];
	fm.template SetIndex<true,false>(updated.index); // useFactoring, smallCorrection

	const NormType norm = GetGuess(updated.index);
	const auto props = norm.Props();
	auto selling = norm.Selling(props.tMin);
	static const int SymDim = SymmetricMatrixType::InternalDimension;
	
	// Prepare the array for holding the neighbor values
	NeighborValuesType values;
	constexpr ScalarType inf = std::numeric_limits<ScalarType>::infinity();
	values.fill(-inf);
	
	// Find the minimizer on each interval
	ScalarType t0 = props.tMin;
	const ScalarType objSign = -props.optimDirection; // Turns maximization into min.
		
	ScalarType updateValue = inf;
	ScalarType tActive = inf; bool isInterior=false;
	while(true){ // Enumerates all optimization intervals
		ScalarType t1 = t0;
		const auto nextStep = selling.NextStep(t1);
		t1=std::min(t1,props.tMax);
				
		// Minimize over interval [t0,t1]
		do{
			// Fill the missing neighbor values
			for(int r=0; r<SymDim; ++r){
				if(values[r]==-inf){
					int orderp=1,orderm=1;
					const auto offset = OffsetType::CastCoordinates(selling.offsets[r]);
					const ScalarType
					valp = fm.template GetNeighborValue<true,false,1>( offset,orderp),
					valm = fm.template GetNeighborValue<true,false,1>(-offset,orderm);
					values[r] = std::min(orderp ? valp : inf, orderm ? valm : inf);
				}
			}
//TODO : optimization opportunity, sort values only once
//TODO : optimization opportunity, get val0 from previous val1
			using Vec1 =  LinearAlgebra::Vector<ScalarType,1>;
			using Diff1 = LinearAlgebra::DifferentiationType<ScalarType,Vec1>;
			using Diff2 = LinearAlgebra::AD2<ScalarType,1>;
			// Note that Diff1(x,0) and Diff2(x,0) differentiate along first coordinate
			
			// Get the update value, and derivative at t0 endpoint
			const Diff1 val0 =
			objSign*norm.UpdateValue(Diff1(t0,0), values, selling);
			if(val0.v[0]>=0 && val0.s<inf){ // Minimum over [t0,t1] attained at t0
				if(val0<updateValue){
					updateValue = val0.s;
					tActive=t0; isInterior=false;
				}
				break;
			}
			
			const Diff1 val1 =
			objSign*norm.UpdateValue(Diff1(t1,0), values, selling);
			if(val1.v[0]<=0 && val1.s<inf){ // Minimum over [t0,t1] attained at t1
				if(val1<updateValue){
					updateValue=val1.s;
					tActive=t1; isInterior=false;
				}
				break;
			}
			
			// There is a degenerate case where the weights associated to all
			// the accepted points vanish (usually a single one). We bypass it here.
			if(std::isinf(val0.s) && std::isinf(val1.s)){continue;}
			
			// Now, root finding in the interval [t0,t1]
			ScalarType x0=t0,x1=t1, // Lower and upper bounds for the subinterval
			x=(val1.v[0]*t0-val0.v[0]*t1)/(val1.v[0]-val0.v[0]),xOld=inf; // Candidate point
			assert(std::isfinite(x));
			while(std::abs(xOld-x)>100*std::numeric_limits<ScalarType>::epsilon()){
				xOld=x;
				// Mathematically guaranteed, but may fail at machine precision
				assert(x0-1e-9<=x && x<=x1+1e-9); // assert(x0<=x && x<=x1);
				const Diff2 val =
				objSign*norm.UpdateValue(Diff2(x,0), values,selling);
				// Update the sub-interval
				if(val.v[0]<=0) {x0=x;}
				else {x1=x;}
				// Try a Newton step
				if(val.m(0,0)>0){
					x-=val.v[0]/val.m(0,0);
					if(x0<x && x<x1) continue;
				}
				// Otherwise do binary subdivision
				x=(x0+x1)/2.;
			}
			assert(std::isfinite(x));
			const ScalarType val = objSign*norm.UpdateValue(x, values,selling);
			if(val<updateValue){
				updateValue=val;
				tActive=x; isInterior=true;
			}
		} while(false);
		
		// Prepare for next interval
		if(t1==props.tMax) break;
		t0 = t1;
		selling.MakeStep(nextStep);
		
		// Update the accepted position
		const auto [i,j] = nextStep;
		const int r = ReductionType::LinearizeIndices(i,j);
		values[r] = -inf;

		// Update the values
		if constexpr(Dimension==3){
			const auto [k,l] = ReductionType::ComplementIndices(i,j);
			std::swap(values[ReductionType::LinearizeIndices(i,k)],
					  values[ReductionType::LinearizeIndices(i,l)]);
		}
	}

	const ScalarType newValue = objSign*updateValue;
	if(newValue<oldValue) {
		active = ActiveNeighFlagType(tActive,isInterior);
		return newValue;
	} else {
		return oldValue;
	}
}

template<int VD> auto StencilTTI<VD>::
_HopfLaxRecompute(IndexCRef index, ActiveNeighFlagType active,
				 DiscreteFlowType & discreteFlow) -> RecomputeType {
	// Possible improvement : introduce some Newton iterations to refine tActive ?
	// (Not strictly necessary to achieve second order, thanks to envelope theorem)
	assert(discreteFlow.empty());
	
	constexpr int SymDim = SymmetricMatrixType::InternalDimension;
	constexpr ScalarType inf = std::numeric_limits<ScalarType>::infinity();
	RecomputeType result;

	// In case of strong anisotropy in a corner.
	if(active.none()) {result.value=inf; result.width=inf; return result;}
	
	auto & fm = *(this->pFM);
	fm.template SetIndex<true,false>(index);
		
	const NormType norm = GetGuess(index);
	const ScalarType t = active.tActive();
	const auto selling = norm.Selling(t);

	NeighborValuesType values;
	std::array<int,SymDim> orders; orders.fill(-1);
	std::array<int,SymDim> signs;
	
	// Get the neighbor values
	for(int r=0; r<SymDim; ++r){
		int orderp=fm.order,orderm=fm.order;
		const OffsetType offset = OffsetType::CastCoordinates(selling.offsets[r]);
		ScalarType
		valp = fm.template GetNeighborValue<true,false,3>( offset,orderp),
		valm = fm.template GetNeighborValue<true,false,3>(-offset,orderm);
		if(orderp==0){valp=inf;}
		if(orderm==0){valm=inf;}
		
		if(valp<valm){
			values[r]=valp;
			orders[r]=orderp;
			signs[r]=1;
		} else {
			values[r]=valm;
			orders[r]=orderm;
			signs[r]=-1;
		}
	}
	constexpr std::array<ScalarType,4> multOrder
	= {0.,1.,3/2.,11/6.};
	 
	// Following code is close to norm.UpdateValue,
	// - but handles high order
	// - does not handle automatic differentiation
	
	// Sort the neighbor values
	std::array<ScalarType,SymDim> indices;
	for(int i=0; i<SymDim; ++i) {indices[i]=i;}
	std::sort(indices.begin(),indices.end(),
			  [&values](int i,int j){return values[i]<values[j];});
	const ScalarType valMin=values[indices[0]];
	assert(valMin<inf);
	
	// Compute the update value
	ScalarType a(0.), b(0.), c(-norm.Multiplier(t)), sol(inf);
	int r=0;
	for(; r<SymDim; ++r){
		const int i=indices[r];
		const ScalarType v=values[i]-valMin;
		if(v>=sol) break;
		const ScalarType w0 = (1-t)*selling.weights0[i]+t*selling.weights1[i];
		const ScalarType w = w0 * square(multOrder[orders[i]]),
		wv=w*v,wvv=wv*v;
		a+=w;
		b+=wv;
		c+=wvv;
		
		if(a <= 100*std::abs(c)*std::numeric_limits<ScalarType>::epsilon()){continue;}
		const ScalarType delta = b*b-a*c;
		assert(delta>=0.);
		sol = (b+sqrt(delta))/a;
	}
	
	// Fill the discrete flow weights
	result.value = sol+valMin;
	result.width=0.;
	ScalarType weightSum = 0.;
	for(int r_=0; r_<r; ++r_){
		const int i=indices[r_];
		const ScalarType v = values[i]-valMin;
		const ScalarType w0 = (1-t)*selling.weights0[i]+t*selling.weights1[i];
		assert(w0>=-1e-10);
		assert(sol-v>=-1e-10);
		weightSum += w0;
		const ScalarType wd = std::max(0.,(sol-v) * w0);
		result.width += wd;
		const ScalarType w = wd * multOrder[orders[i]]; // No square multOrder
		const OffsetType offset = OffsetType::CastCoordinates(selling.offsets[i]);
		discreteFlow.push_back({offset*signs[i],w});
	}
	result.width/=weightSum;
	return result;
}
#endif /* TTI_hpp */
