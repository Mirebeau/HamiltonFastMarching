//
//  Seismic2.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 30/01/2019.
//

#ifndef Seismic2_hpp
#define Seismic2_hpp

template<typename T> auto StencilGenericLag2<T>::
GetNorm(IndexCRef index) const -> NormType {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm((*pMetric)(index), param.gridScale);
}

template<typename T> auto StencilGenericLag2<T>::
GetGuess(const PointType & p) const -> DistanceGuess {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm(MapWeightedSum<MetricElementType>(*pMetric,this->pFM->dom.Neighbors(p)), param.gridScale);
//	const ScalarType invh2 = 1./square(param.gridScale);
//	return NormType{invh2*MapWeightedSum<MetricElementType>(*pMetric,pFM->dom.Neighbors(p))};
}

template<typename T> auto StencilGenericLag2<T>::
HopfLaxUpdate(IndexCRef index, const OffsetVal3 & offsetVal)
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
		value = norm.HopfLax({neigh(0)},Vec<1>{val(0)}).first;
		active = 0;
	}
	
	if(offsetVal.size()>=2){
		value = norm.HopfLax({neigh(0),neigh(1)},Vec<2>{val(0),val(1)}).first;
		active = 1;
	}
	
	if(offsetVal.size()==3){
		const auto hl = norm.HopfLax({neigh(0),neigh(2)},Vec<2>{val(0),val(2)});
		if(hl.first<value){
			value = hl.first;
			active = 2;
		}
	}
	
	return {value,active};
}

template<typename T> auto StencilGenericLag2<T>::
HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow)
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
		const auto & [value,weights] = norm.HopfLax({neigh(0)},Vec<1>{w(0)});
		w(0)=weights[0];
		return {value,0.};
	} else {
		assert(flow.size()==2);
		const auto & [value,weights] = norm.HopfLax({neigh(0),neigh(1)},Vec<2>{w(0),w(1)});
		const ScalarType width = weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1));
		w(0)=weights[0]; w(1)=weights[1];
		assert(weights.Sum()>0);
		return {value,width/weights.Sum()};
	}
}

template<typename T> void StencilGenericLag2<T>::
SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil){
	const NormType & norm = GetNorm(index);
	assert(tmp_stencil.empty());
	tmp_stencil.insert(tmp_stencil.end(),{OffsetType(1,0),OffsetType(0,-1),OffsetType(-1,0),OffsetType(0,1)});
	
	/*
	// Predicate based version, uses two to three gradient evaluations per point.
	auto pred = [&norm,this](OffsetCRef u, OffsetCRef v) -> bool {
		return CosAngle(norm,VectorType::CastCoordinates(u),
						VectorType::CastCoordinates(v)) >= this->cosAngleMin;};
	
	SternBrocotRefine(pred, stencil, tmp_stencil);
	*/
	
	SternBrocotRefine_AcuteBound(norm, cosAngleMin, stencil, tmp_stencil, tmp_stencil_vec, tmp_stencil_scal);
}

template<typename T> void StencilGenericLag2<T>::
Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	pMetric = that->template GetField<MetricElementType>("metric",false);
	cosAngleMin = that->io.Get("cosAngleMin", cosAngleMin);
}
#endif /* Seismic2_h */
