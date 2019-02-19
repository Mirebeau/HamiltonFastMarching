//
//  Seismic2.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 30/01/2019.
//

#ifndef Seismic2_hpp
#define Seismic2_hpp

auto StencilSeismic2::GetNorm(IndexCRef index) const -> NormType {
	const ScalarType invh2 = 1./square(param.gridScale);
	assert(pMetric!=nullptr);
	return NormType{invh2*(*pMetric)(index)};
}

auto StencilSeismic2::HopfLaxUpdate(IndexCRef index, const OffsetVal3 & offsetVal)
-> std::pair<ScalarType,int> {
	assert(!offsetVal.empty());
	const NormType  & norm = GetNorm(index);

	auto neigh = [&offsetVal](int i) -> VectorType {
		assert(i<offsetVal.size());
		return VectorType::CastCoordinates(offsetVal[i].first);
	};
	auto val = [&offsetVal](int i) -> ScalarType {
		assert(i<offsetVal.size());
		return offsetVal[i].second;
	};
	
	int active;
	ScalarType value;

	if(offsetVal.size()==1) {
		value = norm.HopfLax({neigh(0)},{val(0)}).first;
		active = 0;
	}
	
	if(offsetVal.size()>=2){
		value = norm.HopfLax({neigh(0),neigh(1)},{val(0),val(1)}).first;
		active = 1;
	}
	
	if(offsetVal.size()==3){
		const auto hl = norm.HopfLax({neigh(0),neigh(2)},{val(0),val(2)});
		if(hl.first<value){
			value = hl.first;
			active = 2;
		}
	}
	
	return {value,active};
}

auto StencilSeismic2::HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow)
-> RecomputeType {
	assert(!flow.empty());
	const NormType & norm = GetNorm(index);

	auto neigh = [&flow](int i) -> VectorType {
		assert(i<flow.size());
		return VectorType::CastCoordinates(flow[i].offset);
	};
	
	// Initially provides value at neighbor, then stores weight
	auto w = [&flow](int i) -> ScalarType & {
		assert(i<flow.size());
		return flow[i].weight;
	};

	
	if(flow.size()==1){
		const auto & [value,weights] = norm.HopfLax({neigh(0)},{w(0)});
		w(0)=weights[0];
		return {value,0.};
	} else {
		assert(flow.size()==2);
		const auto & [value,weights] = norm.HopfLax({neigh(0),neigh(1)},{w(0),w(1)});
		const ScalarType width = weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1));
		assert((weights.Sum()-1.)<1e-6);
		w(0)=weights[0]; w(1)=weights[1];
		return {value,width};
	}
}

void StencilSeismic2::SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil){
	const NormType & norm = GetNorm(index);
	l.clear();
	l.insert_after(l.before_begin(),{OffsetType(1,0),OffsetType(0,1),OffsetType(-1,0),OffsetType(0,-1),OffsetType(1,0)});
	
	SternBrocotRefine([&norm,this](OffsetCRef u, OffsetCRef v) -> bool {
		return norm.CosAngle(VectorType::CastCoordinates(u), VectorType::CastCoordinates(v)) >= this->cosAngleMin;}, l);

/*	std::cout << index << std::endl;
	for(auto & v : l) std::cout << VectorType::CastCoordinates(v); std::cout << std::endl;*/
	stencil.insert(stencil.end(),l.begin(),l.end());
	stencil.pop_back();
}

void StencilSeismic2::Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	pMetric = that->template GetField<MetricElementType>("metric",false);
	cosAngleMin = that->io.Get("cosAngleMin", cosAngleMin);
}
#endif /* Seismic2_h */
