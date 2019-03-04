//
//  HFM_StencilLag3.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 18/02/2019.
//

#ifndef HFM_StencilLag3_hxx
#define HFM_StencilLag3_hxx

// ********** Stencil data - Semi-Lagrangian3 *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Lag3,Dummy>
: HamiltonFastMarching<T>::_StencilDataTypeBase {
	using HFM = HamiltonFastMarching<T>;
	using Superclass = HFM::_StencilDataTypeBase;
	
	IndexType dims; // Needs value
	virtual void Setup(HFMI *) override;
	virtual void Initialize(const HFM *) override final;
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override final;
	template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
	
	typedef CappedVector<std::pair<OffsetType,ScalarType>, 9> OffsetVals;
	virtual std::pair<ScalarType,int> HopfLaxUpdate(IndexCRef,const OffsetVals &) = 0;
	virtual RecomputeType HopfLaxRecompute(IndexCRef,DiscreteFlowType &) = 0;
private:
	HFM::Array<StencilType,Dimension> stencils;
	std::vector<OffsetType> tmpHLU;
};

// -------------- Implementation --------


// Setup and initialization
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
Setup(HFMI * that){
	Superclass::Setup(that);
	stencils.dims = dims;
}

// Compute the direct and reversed offsets
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
Initialize(const HFM * pFM_){
	Superclass::Initialize(pFM_);
	const DiscreteType size = dims.Product();
	stencils.resize(size);
	std::vector<std::pair<DiscreteType,OffsetType> > targets;
	std::vector<OffsetType> directOffsets;
	
	for(DiscreteType linearIndex=0; linearIndex<size; ++linearIndex){
		const IndexType index = stencils.Convert(linearIndex);
		StencilType & stencil = stencils[linearIndex];
		this->SetStencil(index,stencil);
		
		directOffsets.clear();
		stencil.Neighbors(directOffsets);
		for(const OffsetType & offset : directOffsets){
			IndexType neighbor = index+IndexDiff::CastCoordinates(offset);
			if(this->pFM->dom.Periodize(neighbor,index).IsValid()){
				targets.push_back({stencils.Convert(neighbor), offset});}
		}
	}
	std::sort(targets.begin(),targets.end());
	assert(this->reversedOffsets.empty());
	this->reversedOffsets.insert(targets,this->dims.Product());
}


// HofLaxUpdate
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
HopfLaxUpdate(FullIndexCRef full, OffsetCRef acceptedOffset, ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
	const IndexType updatedIndex = full.index;
	
	// Get the index of the accepted offset
	const StencilType & stencil = stencils[full.linear];
	const ShortType neighborIndex = stencil.NeighborIndex(acceptedOffset);
	
	// Iterate over neighbors, and fill relative data
	OffsetVals offsetVal;
	const HFM & fm = *this->pFM;
	
	tmpHLU.clear();
	stencil.NeighborsAround(neighborIndex, tmpHLU);
	for(auto & offset : tmpHLU){
		IndexType neigh = updatedIndex+IndexDiff::CastCoordinates(offset);
		const auto transform = fm.dom.Periodize(neigh,updatedIndex);
		const ScalarType value =
		( transform.IsValid() && fm.acceptedFlags(neigh) ) ?
		fm.values(neigh) : std::numeric_limits<ScalarType>::infinity();
		offsetVal.push_back({offset,value});
	}
	
	// Put the data of the accepted value
	offsetVal.push_back({acceptedOffset, acceptedValue});
	
	// Evaluate the Hopf-Lax update, and compare
	const auto [value,sectorIndex] = HopfLaxUpdate(updatedIndex,offsetVal);
	const ScalarType oldValue = fm.values[full.linear];
	
	if(value>=oldValue) return oldValue;
	
	active.neighborIndex = neighborIndex;
	active.sectorIndex = sectorIndex;
	return value;
}

// HopfLaxRecompute
template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag3, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef updatedIndex, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
	assert(!active.none());
	const StencilType & stencil = stencils(updatedIndex);
	
	// Get the offsets toward the active neighbors
	// TODO : check if these array allocations are perf limiting, and act if necessary.
	tmpHLU.clear();
	stencil.Neighbors(tmpHLU);
	const OffsetType offset0 = tmpHLU[active.neighborIndex];
	
	tmpHLU.clear();
	stencil.NeighborsAround(active.neighborIndex,tmpHLU);
	const OffsetType offset1 = tmpHLU[active.sectorIndex];
	const OffsetType offset2 = tmpHLU[(active.sectorIndex+1)%tmpHLU.size()];
	
	// Get the values at the active neighbors, and their approximation order
	int ord0 = 3;
	ScalarType val0 = f(offset0,ord0);
	assert(ord0!=0); // This value must be valid, hence with positive approx order
	int ord1 = ord0;
	ScalarType val1 = f(offset1,ord1);
	int ord2 = ord1!=0 ? ord1 : ord0;
	ScalarType val2 = f(offset2,ord2);
	
	// All neighbor values must be at the same approximation order
	const int ord = ord2!=0 ? ord2 : ord1!=0 ? ord1 : ord0;
	if(ord<ord0) {ord0=ord; val0 = f(offset0,ord0); assert(ord0==ord);}
	if(ord<ord1) {ord1=ord; val1 = f(offset1,ord1); assert(ord1==ord);}
	assert(ord==ord2 || ord2==0);
	
	assert(discreteFlow.empty());
	const std::array<ScalarType,4> mult = {0.,1.,3./2.,11./6.};
	discreteFlow.push_back({offset0,mult[ord]*val0});
	if(ord1>0) {discreteFlow.push_back({offset1,mult[ord]*val1});}
	if(ord2>0) {discreteFlow.push_back({offset2,mult[ord]*val2});}
	
	RecomputeType result = HopfLaxRecompute(updatedIndex,discreteFlow);
	const std::array<ScalarType,4> div = {0.,1.,2./3.,6./11.};
	result.width*=div[ord]; result.value*=div[ord];
	
	return result;
};

#endif /* HFM_StencilLag3_hxx */
