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
}

template<int VD> auto StencilTTI<VD>::
MakeNorm(const MetricElementType & m) const -> NormType {
	const auto & [alg,trans] = m;
	const auto & [linear,quadratic] = alg;
	NormType norm{linear,quadratic,trans};
	for(int i=0; i<Dimension; ++i){
		for(int j=0; j<Dimension; ++j){
			norm.transform(i,j)/=param.gridScales[j];
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
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef acceptedOffset, ScalarType value,
			  ActiveNeighFlagType & active) -> ScalarType {
	
	auto & fm = *(this->pFM);
	const ScalarType oldValue = fm.values[updated.linear];
	fm.template SetIndex<true,false>(updated.index); // useFactoring, smallCorrection

	const NormType norm = GetGuess(updated.index);
	const auto props = norm.Props();
	auto selling = norm.Selling(props.tMin);
	
	// Test if the provided offset is within the stencil
	int acceptedPos = -1;
	using DVec = DiscreteVectorType;
	const DVec offsetp = DVec::CastCoordinates(acceptedOffset);
	const DVec offsetm = -offsetp;
	
	static const int SymDim = SymmetricMatrixType::InternalDimension;
	for(int i=0; i<SymDim; ++i){
		const DVec & e = selling.offsets[i];
		if(e==offsetp || e==offsetm){
			acceptedPos = i;
			break;
		}
	}
	
/*	const bool display = (updated.index==IndexType{0,2});
	if(display){
	std::cout << "In HopfLaxUpdate"
	ExportVarArrow(updated.index)
	ExportVarArrow(acceptedOffset)
	ExportVarArrow(value)
	ExportVarArrow(acceptedPos)
	ExportVarArrow(norm) << std::endl;
	}*/
	
	// Prepare the array for holding the neighbor values
	NeighborValuesType values;
	const ScalarType inf = std::numeric_limits<ScalarType>::infinity();
	values.fill(-inf);
	
	// Find the minimizer on each interval
	ScalarType t0 = props.tMin;
	const ScalarType objSign = -props.optimDirection; // Turns maximization into min.
	ScalarType updateValue = inf;
	ScalarType tActive = inf; bool isInterior=false;
	while(true){
		ScalarType t1 = t0;
		const auto nextStep = selling.NextStep(t1);
		t1=std::min(t1,props.tMax);
		
		// Minimize over interval [t0,t1]
		do{
			if(acceptedPos==-1) break;
			// Fill the missing neighbor values
			for(int r=0; r<SymDim; ++r){
				if(values[r]==-inf){
					int orderp=1,orderm=1;
					const auto offset = OffsetType::CastCoordinates(selling.offsets[r]);
					const ScalarType
					valp = fm.template GetNeighborValue<true,false,1>( offset,orderp),
					valm = fm.template GetNeighborValue<true,false,1>(-offset,orderm);
					values[r] = std::min(orderp ? valp : inf, orderm ? valm : inf);
//					assert(std::max(orderp,orderm)==1);
/*					std::cout
					ExportVarArrow(offset)
					ExportVarArrow(valp)
					ExportVarArrow(valm)
					<< std::endl;*/
					
				}
			}
/*			std::cout
			ExportArrayArrow(values) << "\n"
			ExportArrayArrow(fm.values)
			<< std::endl;*/
//TODO : optimization opportunity, sort values only once
//TODO : optimization opportunity, get val0 from previous val1
			using Vec1 =  LinearAlgebra::Vector<ScalarType,1>;
			using Diff1 = LinearAlgebra::DifferentiationType<ScalarType,Vec1>;
			using Diff2 = LinearAlgebra::AD2<ScalarType,1>;
			// Note that Diff1(x,0) and Diff2(x,0) differentiate along first coordinate
			
			// Get the update value, and derivative at t0 endpoint
			const Diff1 val0 =
			objSign*norm.UpdateValue(Diff1(t0,0), values,selling);
//			std::cout ExportVarArrow(val0) ExportVarArrow(t0) << std::endl;
			if(val0.v[0]>=0){ // Minimum over [t0,t1] attained at t0
				if(val0<updateValue){
					updateValue = val0.s;
					tActive=t0; isInterior=false;
				}
				break;
			}
			// Get the update value, and derivative at t0 endpoint
			const Diff1 val1 =
			objSign*norm.UpdateValue(Diff1(t1,0), values, selling);
//			std::cout ExportVarArrow(val1) ExportVarArrow(t1) << std::endl;
			if(val1.v[0]<=0){ // Minimum over [t0,t1] attained at t1
				if(val1<updateValue){
					updateValue=val1.s;
					tActive=t1; isInterior=false;
				}
				break;
			}
			
			{
				std::vector<ScalarType> updates;
				const int nI=100;
				for(int i=0; i<nI; ++i){
					ScalarType x=t0+(t1-t0)*(i/double(nI));
					updates.push_back(norm.UpdateValue(x,values,selling));
				}
//				std::cout ExportArrayArrow(updates) << std::endl;
			}
			
//			std::cout ExportVarArrow(val0) ExportVarArrow(val1) << std::endl;
			// Now, root finding in the interval [t0,t1]
			ScalarType x0=t0,x1=t1, // Lower and upper bounds for the subinterval
			x=(val1.v[0]*t0-val0.v[0]*t1)/(val1.v[0]-val0.v[0]),xOld=inf; // Candidate point
			assert(x0<x && x<x1);
			while(std::abs(xOld-x)>100*std::numeric_limits<ScalarType>::epsilon()){
				xOld=x;
				assert(x0<=x && x<=x1);
				const Diff2 val =
				objSign*norm.UpdateValue(Diff2(x,0), values,selling);
//				std::cout ExportVarArrow(x) ExportVarArrow(val) << std::endl;
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
		if(acceptedPos==r){acceptedPos=-1;}
		const DVec & e = selling.offsets[r];
		if(acceptedPos==-1 && (e==offsetp || e==offsetm)){acceptedPos=r;}
		
		// Update the values
		if constexpr(Dimension==3){
			const auto [k,l] = ReductionType::ComplementIndices(i,j);
			std::swap(values[ReductionType::LinearizeIndices(i,k)],
					  values[ReductionType::LinearizeIndices(i,l)]);
		}
	}
	
/*	std::cout
	ExportVarArrow(updated.index)
	ExportVarArrow(objSign*updateValue)
	<< "\n---------------------\n";*/
	
	const ScalarType newValue = objSign*updateValue;
	if(newValue<oldValue) {
		active = tActive;
		return newValue;
	} else {
		return oldValue;
	}
}

template<int VD> template<typename F> auto StencilTTI<VD>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active,
				 DiscreteFlowType & discreteFlow) -> RecomputeType {
	// TODO : maybe do not exclude, in case of strong anisotropy in a corner ?
	assert(!active.none());
	assert(discreteFlow.empty());
	
	const NormType norm = DistanceGuess(index);
	const auto props = norm.Props();
	const auto selling = norm.Selling(active.tActive);
	
	static const int SymDim = SymmetricMatrixType::InternalDimension;
	NeighborValuesType values;
	int order=3;
	std::array<int,SymDim> orders; orders.fill(-1);
	std::array<int,SymDim> signs;
	
	// Get the neighbor values
	for(int r=0; r<SymDim; ++r){
		int orderp=order,orderm=order;
		ScalarType
		valp = f( selling.offsets[r],orderp),
		valm = f(-selling.offsets[r],orderm);
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
	// TODO : handle high order
	// TODO : introduce some Newton iterations to refine tActive ?
	assert(order==1);
	const ScalarType t = active.tActive();
	const ScalarType val = norm.UpdateValue(t, values,props,selling);
	
	RecomputeType result;
	
	result.value = val;
	result.width = 0.; ScalarType weightSum=0.;
	for(int r=0; r<SymDim; ++r){
		if(val>values[r]){
			const ScalarType weight = selling.weights0[r]*(1.-t)+selling.weights1[r]*t;
			const ScalarType diff = val-values[r];
			discreteFlow.push_back({selling.offsets[r]*signs[r],weight*diff});
			weightSum+=weight; result.diff+=weight*diff;
		}
	}
	result.diff/=weightSum;
	return val;
}
#endif /* TTI_hpp */
