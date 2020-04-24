//
//  HFM_StencilRecompute.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/02/2019.
//

#ifndef HFM_StencilRecompute_hxx
#define HFM_StencilRecompute_hxx

/*
 This class implements a stencil data structure for the Hamilton-Fast-Marching class,
 with the following features :
 - Stencils have a fixed size, and involve weights and offsets. (Eulerian type.)
 - Each point has an independent stencil.
 - They are computed initially, to get the reversed stencils.
 - They are recomputed when the point is first considered, and then cached until it is accepted.
 */

// ******** Stencil data - recompute *********

template<typename T> template<typename Dummy>
struct HamiltonFastMarching<T>::_StencilDataType<SSP::Recomp,Dummy>
: HamiltonFastMarching<T>::_StencilDataTypeBase {
	typedef HamiltonFastMarching<T> HFM;
	using Superclass = HFM::_StencilDataTypeBase;
	Redeclare6Types(HFM,IndexCRef,OffsetCRef,StencilType,DifferenceType,Traits,FullIndexCRef)
	Redeclare5Types(HFM,ParamInterface,HFMI,MultiplierType,DomainType,DistanceGuess)
	Redeclare6Types(Traits,DiscreteType,ScalarType,PointType,VectorType,IndexType,OffsetType)
	Redeclare1Type(StencilType,QuadType)
	Redeclare2Constants(Traits,Dimension,mathPi)
	
	struct RecomputeDataType {StencilType stencil; MultiplierType mult;}; // mult is dummy here
	RecomputeDataType RecomputeData(IndexCRef);
protected:
	friend struct HamiltonFastMarching<Traits>;
	virtual void EraseCache(DiscreteType index) override final {shallowStencilQuads.erase(index);}
	
	virtual ScalarType HopfLaxUpdate(FullIndexCRef, OffsetCRef, ScalarType, ActiveNeighFlagType &) override final;
	template<typename F> RecomputeType HopfLaxRecompute(const F &, IndexCRef, ActiveNeighFlagType, DiscreteFlowType &);	
	virtual void Initialize(const HFM *) override;
private:
	ShallowMap<DiscreteType, std::pair<StencilType, QuadType> > shallowStencilQuads;
};

// ---- HopfLaxUpdate ----

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef offset,
			  ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
	const bool found = shallowStencilQuads.find(updated.linear);
	auto & [stencil,quad] = shallowStencilQuads[updated.linear];
	if(!found) {this->SetStencil(updated.index,stencil);}
	
	return stencil.HopfLaxUpdate(offset,acceptedValue,MultiplierType{},quad,active);
}

// ---- HopfLaxRecompute ------

template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
	StencilType stencil;
	this->SetStencil(index, stencil);	
	return stencil.HopfLaxRecompute(f,MultiplierType(),active,discreteFlow);
}

// --- Recompute data ---

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
RecomputeData(IndexCRef index) -> RecomputeDataType {
	RecomputeDataType result;
	this->SetStencil(index, result.stencil);
	return result;
}

// --- Initialization  ---

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
Initialize(const HFM * pFM) {
	Superclass::Initialize(pFM);
	shallowStencilQuads.resize(this->dims.Product());
	
	typedef std::pair<DiscreteType,OffsetType> IndexOffsetPair;
	std::vector<IndexOffsetPair> offsets;
	offsets.reserve(shallowStencilQuads.size()*HFM::StencilType::nNeigh);
	IndexType updatedIndex;
	auto InsertOffset = [pFM,&offsets,&updatedIndex](OffsetType offset, ScalarType w){
		if(w==0.) return;
		IndexType acceptedIndex;
		const DomainTransformType & transform =
			pFM->VisibleOffset(updatedIndex,offset,acceptedIndex);
		if(transform.IsValid()){
			transform.PullVector(offset);
			offsets.push_back({pFM->values.Convert(acceptedIndex),offset});
		}
	};

	for(DiscreteType linearIndex=0; linearIndex<pFM->values.size(); ++linearIndex){
		updatedIndex = pFM->values.Convert(linearIndex);
		if(HFM::DomainType::periodizeUsesBase && !pFM->dom.PeriodizeNoBase(updatedIndex).IsValid())
			continue;
		StencilType stencil;
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
	this->reversedOffsets.insert(offsets,this->dims.Product());

}

#endif /* HFM_StencilRecompute_hxx */
