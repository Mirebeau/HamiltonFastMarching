//
//  HFM_StencilLag2.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef HFM_StencilLag2_hxx
#define HFM_StencilLag2_hxx

// --------------- Semi-Lagrangian FM-ASR scheme --------------

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::SetStencil(IndexCRef index, StencilType & stencil){
	const DiscreteType linearIndex = indexConverter.Convert(index);
	const DiscreteType start = directOffsetsSplits[linearIndex], end = directOffsetsSplits[linearIndex+1];
	stencil.pOffsets = &directOffsets[start];
	stencil.nOffsets = end-start;
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::ReversedOffsets(FullIndexCRef full)
-> RangeAccessor<OffsetType*> {
	return
	RangeAccessor<OffsetType*>(&reversedOffsets[reversedOffsetsSplits[full.linear]],
							   &reversedOffsets[reversedOffsetsSplits[full.linear+1]]);
}

// Setup and initialization
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::Setup(HFMI * that){
	dims = IndexType::CastCoordinates( that->io.template Get<PointType>("dims") );
	indexConverter.dims = dims;
}

// Compute the direct and reversed offsets
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::Initialize(const HFM * pFM_){
	pFM = pFM_;
	const DiscreteType size = dims.Product();
	
	directOffsetsSplits.reserve(size+1);
	assert(directOffsets.empty() && directOffsetsSplits.empty());
	std::vector<std::pair<DiscreteType,OffsetType> > targets;
	
	for(DiscreteType linearIndex=0; linearIndex<size; ++linearIndex){
		directOffsetsSplits.push_back((DiscreteType)directOffsets.size());
		const IndexType index = indexConverter.Convert(linearIndex);
		SetNeighbors(index,directOffsets);
		
		
		for(auto it = directOffsets.begin()+directOffsetsSplits.back(); it!=directOffsets.end(); ++it){
			const OffsetType & offset = *it;
			IndexType neighbor = index+IndexDiff::CastCoordinates(offset);
			if(pFM->dom.Periodize(neighbor,index).IsValid()){
				targets.push_back({indexConverter.Convert(neighbor), offset});}
		}
	}
	
	directOffsetsSplits.push_back((DiscreteType)directOffsets.size());
	
	assert(reversedOffsets.empty() && reversedOffsetsSplits.empty());
	reversedOffsetsSplits.reserve(size+1);
	reversedOffsets.reserve(targets.size());
	std::sort(targets.begin(),targets.end());
	
	DiscreteType last = -1;
	for(const auto [index,offset] : targets){
		if(last!=index){
			reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());
			last = index;
		}
		reversedOffsets.push_back(offset);
	}
	reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());
}

// HofLaxUpdate
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
HopfLaxUpdate(FullIndexCRef full, OffsetCRef acceptedOffset, ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
	const IndexType updatedIndex = full.index;
	
	// Get the sector
	StencilType stencil;
	SetStencil(updatedIndex, stencil);
	int i = 0;
	for(;stencil.Sector(i,0)!=acceptedOffset; ++i){ // Find the relevant sector
		assert(i<stencil.NSectors());}
	
	OffsetVal3 offsetVal;
	std::array<ShortType, 3> act;
	act[0] = i;
	offsetVal.push_back({acceptedOffset,acceptedValue});
	
	while(true){
		const OffsetType offset = stencil.Sector(i,1);
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		if(!pFM->dom.Periodize(neigh,updatedIndex).IsValid()) break;
		if(!pFM->acceptedFlags(neigh)) break;
		act[1] = i;
		offsetVal.push_back({offset,pFM->values(neigh)});
		break;
	}
	
	while(true){
		const int j = (i==0 ? stencil.NSectors()-1 : i-1);
		const OffsetType offset = stencil.Sector(j,0);
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		if(!pFM->dom.Periodize(neigh,updatedIndex).IsValid()) break;
		if(!pFM->acceptedFlags(neigh)) break;
		act[offsetVal.size()] = j;
		offsetVal.push_back({offset,pFM->values(neigh)});
		break;
	}
	
	const auto [newValue,newActive] = HopfLaxUpdate(updatedIndex,offsetVal);
	const ScalarType oldValue = pFM->values[full.linear];
	
	if(newValue>=oldValue) return oldValue;
	
	active.sectorIndex = act[newActive];
	return newValue;
}

// HopfLaxRecompute
template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
	assert(!active.none());
	StencilType stencil;
	SetStencil(index, stencil);
	
	// Get data of neighbor 0
	int ord0 = 3;
	const OffsetType offset0 = stencil.Sector(active.sectorIndex,0);
	ScalarType val0 = f(offset0,ord0);
	
	
	// Get data of neighbor 1
	int ord1 = ord0==0 ? 3 : ord0;
	const OffsetType offset1 = stencil.Sector(active.sectorIndex,1);
	const ScalarType val1 = f(offset1,ord1);
	
	if(ord1<ord0){
		ord0=ord1;
		val0 = f(offset0,ord0);
		assert(ord0==ord1);
	}
	
	assert(ord0>=0 && ord1>=0);
	assert(ord0+ord1>0);
	const int ord = std::max(ord0,ord1); // snd order used ?
	
	assert(discreteFlow.empty());
	const std::array<ScalarType,4> mult = {0.,1.,3./2.,11./6.};
	if(ord0>0) {discreteFlow.push_back({offset0,mult[ord]*val0});}
	if(ord1>0) {discreteFlow.push_back({offset1,mult[ord]*val1});}
	
	RecomputeType result = HopfLaxRecompute(index,discreteFlow);
	const std::array<ScalarType,4> div = {0.,1.,2./3.,6./11.};
	result.width*=div[ord]; result.value*=div[ord];
	
	return result;
};


#endif /* HFM_StencilLag2_hxx */
