//
//  HFM_StencilShare.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/02/2019.
//

#ifndef HFM_StencilShare_hxx
#define HFM_StencilShare_hxx

/*
This class implements a stencil data structure for the Hamilton-Fast-Marching class,
with the following features :
- Stencils have a fixed size, and involve weights and offsets. (Eulerian type.)
- The weight can be modulated by a (vector) cost function.
- The stencils depend only on some of the point coordinates,
  and are shared between points with those same coordinates.
- They are computed only once.
*/


template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Share,Dummy>
: HamiltonFastMarching<T>::_StencilDataTypeBase {
	using HFM = HamiltonFastMarching<T>;
	using Superclass = HFM::_StencilDataTypeBase;
	Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,MultiplierType,Traits)
	Redeclare5Types(HFM,ParamInterface,HFMI,DomainType,FullIndexCRef,DistanceGuess)
	Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
	Redeclare1Type(StencilType,QuadType)
	Redeclare2Constants(Traits,Dimension,mathPi)
	Redeclare1Type(Superclass, ConstOffsetRange);
	
	typedef typename HFM::template DataSource<MultiplierType> MultSourceType;
	std::unique_ptr<MultSourceType> pMultSource = nullptr; // Needs assignment
	
	struct RecomputeDataType {const StencilType & stencil;MultiplierType mult;};
	RecomputeDataType RecomputeData(IndexCRef); // Needed in some subclasses of HFM...
	virtual ConstOffsetRange ReversedOffsets(FullIndexCRef) const override final;
	
	virtual void Setup(HFMI *) override;
protected:
	friend struct HamiltonFastMarching<Traits>;
	void EraseCache(DiscreteType index) override final {shallowMultQuads.erase(index);}
	
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override final;
	template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);
	
	virtual void Initialize(const HFM *) override;
private:
	typedef typename HFM::template Array<StencilType, T::nStencilDependencies> StencilArrayType;
	StencilArrayType stencils;
	ShallowMap<DiscreteType, std::pair<MultiplierType, QuadType> > shallowMultQuads;
	typedef typename StencilArrayType::IndexType ShortIndexType;
	IndexType IndexFromShortIndex(const ShortIndexType &) const;
	ShortIndexType ShortIndexFromIndex(IndexCRef) const;
};

// ----- Default setup ------

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<StencilStoragePolicy::Share,Dummy>::
Setup(HFMI*that){
	Superclass::Setup(that);
	if(that->io.HasField("speed") || that->io.HasField("speed_times"))
		pMultSource = that->template GetField<MultiplierType>("speed"); // Speed is an external alias for multiplier
	else {
		typedef typename HFMI::template DataSource_Inverse<MultiplierType> SourceInvType;
		pMultSource = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<MultiplierType>("cost") ) );
	}
}

// ---- HopfLaxUpdate ----
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef offset,
			  ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
	const bool found = shallowMultQuads.find(updated.linear);
	auto & [mult,quad] = shallowMultQuads[updated.linear];
	if(!found) { mult = (*pMultSource)(updated.index);}
	const StencilType & stencil = stencils(ShortIndexFromIndex(updated.index));
	
	return stencil.HopfLaxUpdate(offset,acceptedValue,mult,quad,active);
}

// ---- HopfLaxRecompute ------

template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active,
				 DiscreteFlowType & discreteFlow) -> RecomputeType {
	return stencils(ShortIndexFromIndex(index)).
	HopfLaxRecompute(f,(*pMultSource)(index),active,discreteFlow);
};

// --- Recompute data ---

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
RecomputeData(IndexCRef index) -> RecomputeDataType {
	return RecomputeDataType{stencils(ShortIndexFromIndex(index)),(*pMultSource)(index)};
}

// --- Reversed offsets ---
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share,Dummy>::
ReversedOffsets(FullIndexCRef full) const -> ConstOffsetRange {
	return this->reversedOffsets[stencils.Convert(ShortIndexFromIndex(full.index))];
}

// --- Short index conversion ---
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share,Dummy>::
IndexFromShortIndex(const ShortIndexType & shortIndex) const -> IndexType {
	IndexType index;
	for(int i=0; i<Dimension; ++i) {index[i]=this->dims[i]/2;}
	for(int j=0; j<Traits::nStencilDependencies; ++j){
		index[Traits::stencilDependencies[j]] = shortIndex[j];}
	return index;
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share,Dummy>::
ShortIndexFromIndex(const IndexType & index) const -> ShortIndexType {
	ShortIndexType shortIndex;
	for(int i=0; i<Traits::nStencilDependencies; ++i){
		shortIndex[i] = index[Traits::stencilDependencies[i]];}
	return shortIndex;
}

// --- Initialization  ---
template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
Initialize(const HFM * pFM_) {
	Superclass::Initialize(pFM_);

	if(!this->dims.IsPositive())
		ExceptionMacro("StencilDataType initialization error : Incorrect dims.");
	if(!pMultSource)
		ExceptionMacro("StencilDataType initialization error : Unspecified pMultSource.");
	
	shallowMultQuads.resize(this->dims.Product());
	
	for(int i=0; i<Traits::nStencilDependencies; ++i)
		stencils.dims[i] = this->dims[Traits::stencilDependencies[i]];
	stencils.resize(stencils.dims.Product());
	
	typedef std::pair<DiscreteType,OffsetType> IndexOffsetPair;
	std::vector<IndexOffsetPair> offsets;
	offsets.reserve(stencils.size()*HFM::StencilType::nNeigh);
	IndexType updatedIndex;
	auto InsertOffset = [this,&offsets,&updatedIndex](OffsetType offset, ScalarType w){
		if(w==0.) return;
		IndexType acceptedIndex;
		auto transform = this->pFM->VisibleOffset(updatedIndex,offset,acceptedIndex);
		if(transform.IsValid()){
			transform.PullVector(offset);
			offsets.push_back({
				stencils.Convert(ShortIndexFromIndex(acceptedIndex)),offset});
		}
	};
	
	for(DiscreteType linearIndex=0; linearIndex<stencils.size(); ++linearIndex){
		updatedIndex = IndexFromShortIndex(stencils.Convert(linearIndex));
		StencilType & stencil = stencils[linearIndex];
		this->SetStencil(updatedIndex,stencil);
		for(const auto & diffs : stencil.forward)
			for(const auto & diff : diffs)
				InsertOffset( diff.offset, diff.baseWeight);
		for(const auto & diffs : stencil.symmetric)
			for(const auto & diff : diffs){
				InsertOffset( diff.offset, diff.baseWeight);
				InsertOffset(-diff.offset, diff.baseWeight);
			}
	}
	std::sort(offsets.begin(), offsets.end());
	const auto offsetEnd = std::unique(offsets.begin(), offsets.end());
	offsets.resize(offsetEnd-offsets.begin());
	assert(this->reversedOffsets.empty());
	this->reversedOffsets.insert(offsets,stencils.size());
}

#endif /* HFM_StencilShare_hxx */
