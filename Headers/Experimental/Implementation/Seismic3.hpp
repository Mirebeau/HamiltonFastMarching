//
//  Seismic3.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef Seismic3_hpp
#define Seismic3_hpp

template<typename T> auto StencilGenericLag3<T>::
GridScales() const -> const GridScalesType & {
	if constexpr(std::is_same_v<GridScalesType,ScalarType>) {return param.gridScale;}
	else {return param.gridScales;}
}

template<typename T> auto
StencilGenericLag3<T>::
GetNorm(IndexCRef index) const -> NormType {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm((*pMetric)(index),GridScales());
}

template<typename T> auto
StencilGenericLag3<T>::
GetGuess(const PointType & p) const -> DistanceGuess {
	assert(pMetric!=nullptr);
	return Traits::MakeNorm(MapWeightedSum<MetricElementType>(*pMetric,this->pFM->dom.Neighbors(p)),GridScales());
}


template<typename T> auto
StencilGenericLag3<T>::
HopfLaxUpdate(IndexCRef index, const OffsetVals & offsetVal)
-> std::pair<ScalarType,int> {
	/*
	if(index==IndexType{0,4,2}){
		std::cout << "In seismic3 hl "
		ExportVarArrow(index)
		ExportVarArrow(offsetVal.back().first)
		<< std::endl;
		}*/
	
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
	
	if(useHopfLaxCache){ // Variant caching data
		
		const auto & indexConverter = this->pFM->values;
		const DiscreteType linearIndex = indexConverter.Convert(index);
		
		auto vertexHash = [&linearIndex,&offsetVal,this](int i) -> HashKey {
			return this->hash(linearIndex,offsetVal[i].first);};
		auto edgeHash = [&linearIndex,&offsetVal,this](int i,int j) -> std::pair<HashKey,bool> {
			return this->hash(linearIndex,offsetVal[i].first,offsetVal[j].first);};

		
		// Update from accepted value
		int sectorIndex = 0;
		VectorType cache0;
		ScalarType value = norm.HopfLax({acceptedOffset},{acceptedValue},cache0).first;
		vertexCache.insert({vertexHash(nNeigh),cache0});
		
		
		const int max_size = OffsetVals::max_size();
		CappedVector<VectorType,max_size> _vertexCache;
		CappedVector<VectorType,max_size> _edgeCache;
		CappedVector<ScalarType,max_size> _posCache;
		
		_vertexCache.resize(nNeigh);
		_edgeCache.resize(nNeigh);
		_posCache.resize(nNeigh);
		
		// Updates from edges
		for(int i=0; i<nNeigh; ++i){
			if(val(i)==Traits::Infinity()) continue;
			
			// Get cached data

			assert(vertexCache.count(vertexHash(i))==1);
			_vertexCache[i] = vertexCache.find(vertexHash(i))->second;
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i)}, {acceptedValue,val(i)},
										 {cache0,_vertexCache[i]},_edgeCache[i]);
			
			// Save cached data
			_posCache[i] = hl.second[1]/hl.second.Sum();
			const auto [edgeHashKey,orient] = edgeHash(nNeigh,i);
			edgeCache.insert({edgeHashKey,
				{_edgeCache[i], orient ? _posCache[i] : 1-_posCache[i]}});
							 
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
			
		
		// Updates from faces
		for(int i=0; i<nNeigh; ++i){
			const int j = (i+1)%nNeigh;
			if(val(i)==Traits::Infinity() || val(j)==Traits::Infinity()) continue;
			
			// Get cached data
			const auto [edgeHashKey,orient] = edgeHash(i,j);
			assert(edgeCache.count(edgeHashKey)==1);
			const auto [g_ij,pos_ij] = edgeCache.find(edgeHashKey)->second;
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)},
										 {acceptedValue,val(i),val(j)},
										 {cache0,_vertexCache[i],_vertexCache[j]},
										 {_edgeCache[i],g_ij,_edgeCache[j]},
										 {_posCache[i], (orient ? pos_ij : (1-pos_ij)), 1-_posCache[j]});
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		return {value,sectorIndex};
		
		
	} if(true) { // Recomputation based variant, with a bit of purely local caching.
	
#ifndef XSIMD_HPP
		// Update from accepted value
		int sectorIndex = 0;
		VectorType cache0;
		ScalarType value = norm.HopfLax({acceptedOffset},{acceptedValue},cache0).first;
		
		
		const int max_size = OffsetVals::max_size();
		std::array<VectorType,max_size> _vertexCache;
		std::array<VectorType,max_size> _edgeCache;
		std::array<ScalarType,max_size> _posCache;
		
		// Updates from edges
		for(int i=0; i<nNeigh; ++i){
			if(val(i)==Traits::Infinity()) continue;
			
			// Recompute update for vertex i to set cache
			norm.HopfLax({offset(i)},{val(i)},_vertexCache[i]);
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i)}, {acceptedValue,val(i)},
										 {cache0,_vertexCache[i]},_edgeCache[i]);
			// Save edge cache
			_posCache[i] = hl.second[1]/hl.second.Sum();
			
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		
		// Updates from faces
		for(int i=0; i<nNeigh; ++i){
			const int j = (i+1)%nNeigh;
			if(val(i)==Traits::Infinity() || val(j)==Traits::Infinity()) continue;
			
			// Recompute update for edge (i,j) to set cache
			VectorType g_ij;
			const auto hl_ij = norm.HopfLax({offset(i),offset(j)}, {val(i),val(j)},
										 {_vertexCache[i],_vertexCache[j]},g_ij);
			const ScalarType pos_ij = hl_ij.second[1]/hl_ij.second.Sum();
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)},
										 {acceptedValue,val(i),val(j)},
										 {cache0,_vertexCache[i],_vertexCache[j]},
										 {_edgeCache[i],g_ij,_edgeCache[j]},
										 {_posCache[i], pos_ij, 1-_posCache[j]});
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		return {value,sectorIndex};
#else
		const int max_size = OffsetVals::max_size();
		std::array<VectorType,max_size> _vertexCache;
		std::array<VectorType,max_size> _edgeCache;
		std::array<ScalarType,max_size> _posCache;

		using SimdNormType = typename Traits::SimdNormType;
		using SimdVectorType = typename SimdNormType::VectorType;
		using SimdScalarType = typename SimdNormType::ComponentType;
		constexpr const int simd_size = SimdScalarType::size;

		SimdNormType simd_norm;
		for(int i=0; i<norm.hookeTensor.data.size(); ++i){
			simd_norm.hookeTensor.data[i] = SimdScalarType(norm.hookeTensor.data[i]);}

		SimdVectorType simd_vector;
		CappedVector<int,simd_size> simd_ind;
		
		for(int i=0; i<offsetVal.size(); ++i){
			if(val(i)==Traits::Infinity()) continue;
			const int n = simd_ind.size();
			simd_ind.push_back(i);
			const VectorType offseti = offset(i);
			for(int j=0; j<Dimension; ++j){simd_vector[j][n]=offseti[j];}
			if(simd_ind.size()==simd_size || i==nNeigh){
				if (simd_ind.size()==0) {
				} else if(simd_ind.size()==1){
					_vertexCache[simd_ind[0]] = norm.GradNorm(offset(simd_ind[0]));
				} else {
					simd_vector = simd_norm.GradNorm(simd_vector);
					for(int k=0; k<simd_ind.size(); ++k){
						for(int j=0; j<Dimension; ++j) {
							_vertexCache[simd_ind[k]][j]=simd_vector[j][k];}
					}
					simd_ind.clear();
				}
			}
		}
		
		// Diagostic of optimization opportunities
		/*{
			static int n_vert_simd = 0;
			static int n_vert_tot = 0;
			static int n_edge_simd = 0;
			static int n_edge_tot = 0;
			
			int n_vert=1, n_edge=0;
			for(int i=0; i<offsetVal.size(); ++i){
				if(val(i)==Traits::Infinity()) continue;
				++n_vert;
				++n_edge;
				if(val((i+1)%nNeigh)==Traits::Infinity()) continue;
				++n_edge;
			}
			
			n_vert_simd += (n_vert+3)/4;
			n_vert_tot += n_vert;
			n_edge_simd += (n_edge+3)/4;
			n_edge_tot += n_edge;
			
			static int oldCounter=10;
			static int counter=0;
			++counter;
			if(counter==2*oldCounter){
				oldCounter=counter;
				std::cout << "Simd reduction ratio "
				ExportVarArrow(n_vert_simd/double(n_vert_tot))
				ExportVarArrow(n_edge_simd/double(n_edge_tot))
				<< std::endl;
				}
		}*/
			
		// Update from accepted value
		int sectorIndex = 0;
		const VectorType & cache0 = _vertexCache[nNeigh];
//		ScalarType value = norm.HopfLax({acceptedOffset},{acceptedValue},cache0).first;
		ScalarType value = acceptedValue+acceptedOffset.ScalarProduct(_vertexCache[nNeigh]);
		
		
		// Updates from edges
		for(int i=0; i<nNeigh; ++i){
			if(val(i)==Traits::Infinity()) continue;
			
			// Recompute update for vertex i to set cache
//			norm.HopfLax({offset(i)},{val(i)},_vertexCache[i]);
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i)}, {acceptedValue,val(i)},
										 {cache0,_vertexCache[i]},_edgeCache[i]);
			// Save edge cache
			_posCache[i] = hl.second[1]/hl.second.Sum();
			
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		
		// Updates from faces
		for(int i=0; i<nNeigh; ++i){
			const int j = (i+1)%nNeigh;
			if(val(i)==Traits::Infinity() || val(j)==Traits::Infinity()) continue;
			
			// Recompute update for edge (i,j) to set cache
			VectorType g_ij;
			const auto hl_ij = norm.HopfLax({offset(i),offset(j)}, {val(i),val(j)},
											{_vertexCache[i],_vertexCache[j]},g_ij);
			const ScalarType pos_ij = hl_ij.second[1]/hl_ij.second.Sum();
			
			const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)},
										 {acceptedValue,val(i),val(j)},
										 {cache0,_vertexCache[i],_vertexCache[j]},
										 {_edgeCache[i],g_ij,_edgeCache[j]},
										 {_posCache[i], pos_ij, 1-_posCache[j]});
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		return {value,sectorIndex};
		
#endif
	
	} else { // Never enabled // Straightforward recomputation based variant (costly)
		
		// Update from accepted offset
		const auto hl = norm.HopfLax({acceptedOffset},Vec<1>{acceptedValue});
		ScalarType value = hl.first;
		int sectorIndex = 0;
		
		// Updates from edges
		for(int i=0; i<nNeigh; ++i){
			if(val(i)==Traits::Infinity()) continue;
			const auto hl = norm.HopfLax({acceptedOffset,offset(i)}, Vec<2>{acceptedValue,val(i)});
			
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		// Updates from faces
		for(int i=0; i<nNeigh; ++i){
			const int j = (i+1)%nNeigh;
			if(val(i)==Traits::Infinity() || val(j)==Traits::Infinity()) continue;
			const auto hl = norm.HopfLax({acceptedOffset,offset(i),offset(j)},
										 Vec<3>{acceptedValue,val(i),val(j)});
			if(hl.first>=value) continue;
			value = hl.first;
			sectorIndex = i;
		}
		
		return {value,sectorIndex};
	}

}

template<typename T> auto
StencilGenericLag3<T>::
HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow)
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
		const auto & [value,weights] = norm.HopfLax({offset(0)},Vec<1>{w(0)});
		w(0)=weights[0];
		return {value,0.};
	} else if(flow.size()==2){
		const auto & [value,weights] = norm.HopfLax({offset(0), offset(1)},Vec<2>{w(0),w(1)});
		assert(weights.Sum()>0);
		const ScalarType width = (weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1)))/weights.Sum();
		w(0)=weights[0]; w(1)=weights[1];
		return {value,width};
	} else {
		assert(flow.size()==3);
		const auto & [value,weights] = norm.HopfLax({offset(0), offset(1), offset(2)},Vec<3>{w(0),w(1),w(2)});
		assert(weights.Sum()>0);
		const ScalarType width =
		(weights[0]*abs(value-w(0))+weights[1]*abs(value-w(1))+weights[2]*abs(value-w(2)))/weights.Sum();
		w(0)=weights[0]; w(1)=weights[1]; w(2)=weights[2];
		return {value,width};
	}
	
}


template<typename T> void
StencilGenericLag3<T>::
SetStencil(IndexCRef index, StencilType & stencil){
	// We'll put metric dependent adaptive stencils here in time
	assert(false);
	ExceptionMacro("Seismic3 error : no stencil geometry specified");
	assert(!checkAcuteness);
}

template<typename T> void
StencilGenericLag3<T>::
Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	auto & io=that->io;
	pMetric = that->template GetField<MetricElementType>("metric",false);
	checkAcuteness = (bool)io.template Get<ScalarType>("checkAcuteness",checkAcuteness);
}

// ----- Cache management ------

template<typename T> auto
StencilGenericLag3<T>::
hash(DiscreteType index,OffsetType offset)
-> HashKey {
	HashKey result=index;
	result = (result<<1) +1;
	for(int i=0; i<Dimension; ++i){
		result = (result<<8) + offset[i];}
	return result;
}


template<typename T> auto
StencilGenericLag3<T>::
hash(DiscreteType index,OffsetType offset1, OffsetType offset2)
->  std::pair<HashKey, bool> {
	using Comp = typename OffsetType::LexicographicCompare;
	bool ordered = Comp()(offset1,offset2);
	if(!ordered) std::swap(offset1,offset2);
	
	HashKey result=index;
	result = (result<<1) +1;
	static_assert(sizeof(HashKey)>=8);
	for(int i=0; i<Dimension; ++i){ result = (result<<5) + offset1[i];}
	for(int i=0; i<Dimension; ++i){ result = (result<<5) + offset2[i];}
	
	return {result,ordered};
}

template<typename T> void
StencilGenericLag3<T>::
EraseCache(DiscreteType index) {
	if(!useHopfLaxCache) return;
	
	{ // Erase vertex data
		const HashKey
		lbound = (HashKey(index)<<1) << (8*Dimension),
		ubound = ((HashKey(index)<<1)+2) << (8*Dimension);
		
		/*
		if(lbound<= 8539603201 && 8539603201 <= ubound){
			std::cout << "Erased "
			ExportVarArrow(index)
			<< std::endl;
		}*/
		
		auto & hlCache = vertexCache;
		hlCache.erase(hlCache.lower_bound(lbound),hlCache.lower_bound(ubound));
	}
	
	{ // Erase edge data
		const HashKey
		lbound = (HashKey(index)<<1) << (5*2*Dimension),
		ubound = ((HashKey(index)<<1)+2) << (5*2*Dimension);
		
		auto & hlCache = edgeCache;
		hlCache.erase(hlCache.lower_bound(lbound),hlCache.lower_bound(ubound));
	}
}


#endif /* Seismic3_hpp */
