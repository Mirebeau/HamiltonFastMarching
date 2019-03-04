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
	auto offset = [nNeigh,&offsetVal](int i) -> VectorType {
		assert(i<nNeigh);
		return VectorType::CastCoordinates(offsetVal[i].first);
	};
	auto val = [nNeigh,&offsetVal](int i) -> ScalarType {
		assert(i<nNeigh);
		return offsetVal[i].second;
	};
	
	// Evaluate the operator
	int sectorIndex = -1;
	ScalarType value = std::numeric_limits<ScalarType>::infinity();
	
	// TODO : optimizations, to avoid recomputing edge and vertex related data.
	
	for(int i=0; i<nNeigh; ++i){
		const int j = (i+1)%nNeigh;
		const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)}, {acceptedValue,val(i),val(j)});
		if(hl.first<value){
			value = hl.first;
			sectorIndex = i;
		}
	}
	
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
		const ScalarType width = weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1));
		assert(std::abs(weights.Sum()-1) < 1e-6);
		w(0)=weights[0]; w(1)=weights[1];
		return {value,width};
	} else {
		assert(flow.size()==3);
		const auto & [value,weights] = norm.HopfLax({offset(0), offset(1), offset(2)},{w(0),w(1),w(2)});
		const ScalarType width = weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1))+weights[2]*abs(value-w(2));
		assert(std::abs(weights.Sum()-1) < 1e-6);
		w(0)=weights[0]; w(1)=weights[1]; w(2)=weights[2];
		return {value,width};
	}
	
}


void StencilSeismic3::SetStencil(IndexCRef index, StencilType & stencil){
	// For now, we use the basic cubic stencil.
	// We'll check if it is acute for the different media, in particular argileous ones
	stencil.offsets = { OffsetType(1,0,0), OffsetType(0,1,0), OffsetType(0,0,1) };
	stencil.geom = Lagrangian3StencilGeometry::Cube;
	assert(!checkAcuteness);
}

void StencilSeismic3::Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	pMetric = that->template GetField<MetricElementType>("metric",false);
	checkAcuteness = (bool)that->io.Get<ScalarType>("checkAcuteness",checkAcuteness);
}

#endif /* Seismic3_hpp */
