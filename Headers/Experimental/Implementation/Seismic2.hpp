//
//  Seismic2.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 30/01/2019.
//

#ifndef Seismic2_hpp
#define Seismic2_hpp

template<typename T> auto StencilGenericLag2<T>::
GridScales() const -> const GridScalesType & {
	if constexpr(std::is_same_v<GridScalesType,ScalarType>) {return param.gridScale;}
	else {return param.gridScales;}
}

template<typename T> auto StencilGenericLag2<T>::
GetNorm(IndexCRef index) const -> NormType {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm((*pMetric)(index),GridScales());
}

template<typename T> auto StencilGenericLag2<T>::
GetGuess(const PointType & p) const -> DistanceGuess {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm(
		MapWeightedSum<MetricElementType>(*pMetric,this->pFM->dom.Neighbors(p)),
		GridScales());
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
	
	if constexpr(Traits::useCache==GradientCachingStrategy::Full) {
		// --- Version fully caching the gradients (counter productive it seems) ---
		const DiscreteType linearIndex = this->indexConverter.Convert(index);
		
		auto hash = [&offsetVal,&linearIndex,this](int i) -> long {
			assert(i<offsetVal.size());
			return this->hash(linearIndex,offsetVal[i].first);
		};
		
		// Compute and save gradient at new vertex
		const VectorType neigh0 = neigh(0);
		const ScalarType val0 = val(0);

		// Test from accepted vertex
		VectorType cache0;
		value = norm.HopfLax({neigh0},Vec<1>{val0},cache0).first;
		active = 0;

		// Cache computed data
		auto & vCache = vertexCache;
		const long hash0 = hash(0);
		assert(vCache.count(hash0)==0);
		[[maybe_unused]] const auto insertion_result = vCache.insert({hash0,cache0});
		assert(insertion_result.second);
		vertexCacheKeys.insert({linearIndex,hash0});
		
		VectorType dummyCache;
		
		if(offsetVal.size()>=2){
			// Neighbor has been accepted before, relevant data is cached
			const auto cache1_it = vCache.find(hash(1));
			assert(cache1_it != vCache.end());
			
			value = norm.HopfLax({neigh0,neigh(1)},Vec<2>{val0,val(1)},
			{cache0,cache1_it->second},dummyCache).first;
			active = 1;
		}
		
		if(offsetVal.size()==3){
			// Neighbor has been accepted before
			const auto cache2_it = vCache.find(hash(2));
			assert(cache2_it != vCache.end());
			
			const ScalarType newValue =
			norm.HopfLax({neigh0,neigh(2)},Vec<2>{val0,val(2)},
						 {cache0,cache2_it->second},dummyCache).first;
			if(newValue<value){
				value=newValue;
				active = 2;
			}
		}
		
	} else if constexpr(Traits::useCache == GradientCachingStrategy::Local) {
		// --- Recomputation based variant with a bit of local caching ---
		
		// Compute and save gradient at new vertex
		const VectorType neigh0 = neigh(0);
		const ScalarType val0 = val(0);
		
		// Test from accepted vertex
		VectorType cache0;
		value = norm.HopfLax({neigh0},Vec<1>{val0},cache0).first;
		active = 0;
		
		// Cache computed data
		VectorType dummyCache;
		
		if(offsetVal.size()>=2){
			// Recompute from neighbor, to set cache
			VectorType cache1;
			norm.HopfLax({neigh(1)},Vec<1>{val(1)},cache1);
			
			value = norm.HopfLax({neigh0,neigh(1)},Vec<2>{val0,val(1)},
								 {cache0,cache1},dummyCache).first;
			active = 1;
		}
		
		if(offsetVal.size()==3){
			// Recompute from neighbor, to set cache
			VectorType cache2;
			norm.HopfLax({neigh(2)},Vec<1>{val(2)},cache2);

			const ScalarType newValue =
			norm.HopfLax({neigh0,neigh(2)},Vec<2>{val0,val(2)},
						 {cache0,cache2},dummyCache).first;
			if(newValue<value){
				value=newValue;
				active = 2;
			}
		}
		
	} else if constexpr(Traits::useCache==GradientCachingStrategy::None){
		// --- Version with full gradient recomputation ---
		// Most costly for seismic2 norms, ok if gradient is cheap
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
	// Predicate based version, uses two or three gradient evaluations per point, instead of one.
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
	if(Traits::useCache!=GradientCachingStrategy::Full){
		const IndexType dims = this->indexConverter.dims;
		const size_t front_size_estim = dims.Product()/(*std::max_element(dims.begin(),dims.end()));
		const size_t active_neighbors_estim = 5*front_size_estim;
		vertexCacheKeys.reserve(active_neighbors_estim);
		vertexCache.reserve(active_neighbors_estim);
	}
}

// ----- Cache data ----

template<typename T> long
StencilGenericLag2<T>::
hash(DiscreteType index,OffsetType offset){
	long result=index;
	result = (result<<1) +1;
	for(int i=0; i<Dimension; ++i){
		result = (result<<8) + offset[i];}
	return result;
}

template<typename T> void
StencilGenericLag2<T>::
EraseCache(DiscreteType index) {
	if(Traits::useCache!=GradientCachingStrategy::Full) return;
	const auto rg = vertexCacheKeys.equal_range(index);
	for(auto [ind,key] : RangeAccessor{rg.first,rg.second}) {vertexCache.erase(key);}
	vertexCacheKeys.erase(rg.first,rg.second);
	
	/*
	static int iter = 0;
	++iter;
	const int size = this->indexConverter.dims.Product();
	if(iter%(size/10) == 0) {
		std::cout  ExportVarArrow(vertexCache.size()/double(size)) << std::endl;}
	 */
		/* // Slow map based version
	const long
	lbound = (index<<1) << (8*Dimension),
	ubound = ((index<<1)+2) << (8*Dimension);
	
	auto & hlCache = hopfLaxCache;
	hlCache.erase(hlCache.lower_bound(lbound),hlCache.lower_bound(ubound));
		 */
}
#endif /* Seismic2_h */
