//
//  Seismic3.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef Seismic3_hpp
#define Seismic3_hpp


auto StencilSeismic3::GetNorm(IndexCRef index) const -> NormType {
	const ScalarType invh2 = 1./square(param.gridScale);
	assert(pMetric!=nullptr);
	return NormType{invh2*(*pMetric)(index)};
}

auto StencilSeismic3::GetGuess(const PointType & p) const -> NormType {
	const ScalarType invh2 = 1./square(param.gridScale);
	assert(pMetric!=nullptr);
	return NormType{invh2*MapWeightedSum<MetricElementType>(*pMetric,pFM->dom.Neighbors(p))};
}



auto StencilSeismic3::HopfLaxUpdate(IndexCRef index, const OffsetVals & offsetVal)
-> std::pair<ScalarType,int> {
	assert(!offsetVal.empty());
	const NormType & norm = GetNorm(index);
	
	// Get data of centered offset
	const VectorType acceptedOffset = VectorType::CastCoordinates(offsetVal.back().first);
	const ScalarType acceptedValue = offsetVal.back().second;
	
	// Make accessors to other offsets data
	const int nNeigh = offsetVal.size()-1;
	auto offset = [&offsetVal](int i) -> VectorType {
		return VectorType::CastCoordinates(offsetVal[i].first);};
	auto val = [&offsetVal](int i) -> ScalarType {
		return offsetVal[i].second;};
	
	// Update from accepted value
	const auto hl = norm.HopfLax({acceptedOffset},{acceptedValue});
	ScalarType value = hl.first;
	int sectorIndex = 0;

	// Updates from edges
	for(int i=0; i<nNeigh; ++i){
		if(val(i)==Traits::Infinity()) continue;
		const auto hl = norm.HopfLax({acceptedOffset,offset(i)}, {acceptedValue,val(i)});
		if(hl.first>=value) continue;
		value = hl.first;
		sectorIndex = i;
	}
	
	// Updates from faces
	for(int i=0; i<nNeigh; ++i){
		const int j = (i+1)%nNeigh;
		if(val(i)==Traits::Infinity() || val(j)==Traits::Infinity()) continue;
		const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)}, {acceptedValue,val(i),val(j)});
		if(hl.first>=value) continue;
		value = hl.first;
		sectorIndex = i;
	}
	
	// ! TODO ! : optimizations, to avoid recomputing edge and vertex related data.

	return {value,sectorIndex};
}

auto StencilSeismic3::HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow)
-> RecomputeType {
	assert(!flow.empty());
	const NormType & norm = GetNorm(index);
	
	auto offset = [&flow](int i) -> VectorType {
		assert(i<flow.size());
		return VectorType::CastCoordinates(flow[i].offset);
	};
	
	// Initially provides value at neighbor, then stores weight
	auto w = [&flow](int i) -> ScalarType & {
		assert(i<flow.size());
		return flow[i].weight;
	};
	
	if(flow.size()==1){
		const auto & [value,weights] = norm.HopfLax({offset(0)},{w(0)});
		w(0)=weights[0];
		return {value,0.};
	} else if(flow.size()==2){
		const auto & [value,weights] = norm.HopfLax({offset(0), offset(1)},{w(0),w(1)});
		assert(weights.Sum()>0);
		const ScalarType width = (weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1)))/weights.Sum();
		w(0)=weights[0]; w(1)=weights[1];
		return {value,width};
	} else {
		assert(flow.size()==3);
		const auto & [value,weights] = norm.HopfLax({offset(0), offset(1), offset(2)},{w(0),w(1),w(2)});
		assert(weights.Sum()>0);
		const ScalarType width =
		(weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1))+weights[2]*abs(value-w(2)))/weights.Sum();
		w(0)=weights[0]; w(1)=weights[1]; w(2)=weights[2];
		return {value,width};
	}
	
}


void StencilSeismic3::SetStencil(IndexCRef index, StencilType & stencil){
	// We'll put metric dependent adaptive stencils here in time
	assert(false);
	assert(!checkAcuteness);
}

void StencilSeismic3::Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	auto & io=that->io;
	pMetric = that->template GetField<MetricElementType>("metric",false);
	checkAcuteness = (bool)io.Get<ScalarType>("checkAcuteness",checkAcuteness);
}

#endif /* Seismic3_hpp */
