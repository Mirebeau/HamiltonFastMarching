//
//  HFM_StencilLag2.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef HFM_StencilLag2_hxx
#define HFM_StencilLag2_hxx

// ********** Stencil data - Semi-Lagrangian2 *******

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag2,Dummy> :
HamiltonFastMarching<T>::_StencilDataTypeBase {
	using HFM=HamiltonFastMarching<T>;
	using Superclass = HFM::_StencilDataTypeBase;
	Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,Traits,FullIndexCRef)
	Redeclare5Types(HFM,ParamInterface,HFMI,MultiplierType,DomainType,DistanceGuess)
	Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
	Redeclare1Constant(Traits,Dimension)
	
	void SetStencil(IndexCRef, StencilType &);
	virtual void Setup(HFMI *) override;
	virtual void Initialize(const HFM *) override final;
	
	virtual void SetNeighbors(IndexCRef, std::vector<OffsetType> &) = 0;
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override final;
	template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
	
	typedef CappedVector<std::pair<OffsetType,ScalarType>, 3> OffsetVal3;
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,const OffsetVal3 &) = 0;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) = 0;
private:
	HFM::Array<ScalarType,Dimension> indexConverter;
	MultiVector<OffsetType,DiscreteType> directOffsets;
};

// --------------- Semi-Lagrangian FM-ASR scheme --------------

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::SetStencil(IndexCRef index, StencilType & stencil){
	const auto rg = directOffsets[indexConverter.Convert(index)];
	stencil.pOffsets = &rg.front();
	stencil.nOffsets = rg.end();
}

// Setup and initialization
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
Setup(HFMI * that){
	Superclass::Setup(that);
	indexConverter.dims = this->dims;
}

// Compute the direct and reversed offsets
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
Initialize(const HFM * _pFM){
	Superclass::Initialize(_pFM);
	const DiscreteType size = this->dims.Product();
	
	assert(directOffsets.empty());
	directOffsets.splits.reserve(size+1);
	std::vector<std::pair<DiscreteType,OffsetType> > targets;
	
	for(DiscreteType linearIndex=0; linearIndex<size; ++linearIndex){
		const IndexType index = indexConverter.Convert(linearIndex);
		SetNeighbors(index,directOffsets.values);
		directOffsets.splits.push_back((DiscreteType)directOffsets.values.size());
		
		for(const OffsetType & offset : directOffsets.back()){
			IndexType neighbor = index+IndexDiff::CastCoordinates(offset);
			if(this->pFM->dom.Periodize(neighbor,index).IsValid()){
				targets.push_back({indexConverter.Convert(neighbor), offset});}
		}
	}
	
	std::sort(targets.begin(),targets.end());
	assert(this->reversedOffsets.empty());
	this->reversedOffsets.insert(targets,size);
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
	
	const auto & fm = *this->pFM;
	
	while(true){
		const OffsetType offset = stencil.Sector(i,1);
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		if(fm.dom.Periodize(neigh,updatedIndex).IsValid()) break;
		if(fm.acceptedFlags(neigh)) break;
		act[1] = i;
		offsetVal.push_back({offset,fm.values(neigh)});
		break;
	}
	
	while(true){
		const int j = (i==0 ? stencil.NSectors()-1 : i-1);
		const OffsetType offset = stencil.Sector(j,0);
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		if(fm.dom.Periodize(neigh,updatedIndex).IsValid()) break;
		if(fm.acceptedFlags(neigh)) break;
		act[offsetVal.size()] = j;
		offsetVal.push_back({offset,fm.values(neigh)});
		break;
	}
	
	const auto [newValue,newActive] = HopfLaxUpdate(updatedIndex,offsetVal);
	const ScalarType oldValue = fm.values[full.linear];
	
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
