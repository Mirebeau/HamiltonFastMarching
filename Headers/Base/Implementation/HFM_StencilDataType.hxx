// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef StencilDataType_h
#define StencilDataType_h

// Note : walls must not be taken into account at stencil initialization with mult based.


// ----- Default setup ------

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<StencilStoragePolicy::Share,Dummy>::Setup(HFMI*that){
    dims = IndexType::CastCoordinates( that->io.template Get<PointType>("dims") );
    if(that->io.HasField("speed") || that->io.HasField("speed_times"))
        pMultSource = that->template GetField<MultiplierType>("speed"); // Speed is an external alias for multiplier
    else {
        typedef typename HFMI::template DataSource_Inverse<MultiplierType> SourceInvType;
        pMultSource = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<MultiplierType>("cost") ) );
    }
}

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp,Dummy>::Setup(HFMI*that){
    dims = IndexType::CastCoordinates( that->io.template Get<PointType>("dims") );
}

// ---- HopfLaxUpdate ----
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef offset,
              ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
    const bool found = shallowMultQuads.find(updated.linear);
    auto & multQuad = shallowMultQuads[updated.linear];
    if(!found) { multQuad.first = (*pMultSource)(updated.index);}
    const StencilType & stencil = stencils(ShortIndexFromIndex(updated.index));
    
    return stencil.HopfLaxUpdate(offset,acceptedValue,multQuad.first,multQuad.second,active);
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
HopfLaxUpdate(FullIndexCRef updated, OffsetCRef offset,
              ScalarType acceptedValue, ActiveNeighFlagType & active) -> ScalarType {
    const bool found = shallowStencilQuads.find(updated.linear);
    auto & stencilQuad = shallowStencilQuads[updated.linear];
    if(!found) {SetStencil(updated.index,stencilQuad.first);}
    
    return stencilQuad.first.HopfLaxUpdate(offset,acceptedValue,MultiplierType{},stencilQuad.second,active);
}

// ---- HopfLaxRecompute ------

template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
    return stencils(ShortIndexFromIndex(index)).HopfLaxRecompute(f,(*pMultSource)(index),active,discreteFlow);
};

template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
    StencilType stencil;
    SetStencil(index, stencil);
    return stencil.HopfLaxRecompute(f,MultiplierType(),active,discreteFlow);
}

// --- Recompute data ---

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share, Dummy>::
RecomputeData(IndexCRef index) -> RecomputeDataType {
    return RecomputeDataType{stencils(ShortIndexFromIndex(index)),(*pMultSource)(index)};
}


template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
RecomputeData(IndexCRef index) -> RecomputeDataType {
    RecomputeDataType result;
    SetStencil(index, result.stencil);
    return result;
}

// --- Reversed offsets ---
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share,Dummy>::
ReversedOffsets(FullIndexCRef full)->RangeAccessor<OffsetType*>{
    const DiscreteType i = stencils.Convert(ShortIndexFromIndex(full.index));
    return
    RangeAccessor<OffsetType*>(&reversedOffsets[reversedOffsetsSplits[i]],
                               &reversedOffsets[reversedOffsetsSplits[i+1]]);
}

template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp,Dummy>::
ReversedOffsets(FullIndexCRef full)->RangeAccessor<OffsetType*>{
    return
    RangeAccessor<OffsetType*>(&reversedOffsets[reversedOffsetsSplits[full.linear]],
                               &reversedOffsets[reversedOffsetsSplits[full.linear+1]]);
}

// --- Short index conversion ---
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Share,Dummy>::
IndexFromShortIndex(const ShortIndexType & shortIndex) const -> IndexType {
    IndexType index;
    for(int i=0; i<Dimension; ++i) {index[i]=dims[i]/2;}
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
Initialize(const HFM * pFM) {
    if(!dims.IsPositive())
        ExceptionMacro("StencilDataType initialization error : Incorrect dims.");
    if(!pMultSource)
        ExceptionMacro("StencilDataType initialization error : Unspecified pMultSource.");
        
    shallowMultQuads.resize(dims.Product());
    
    for(int i=0; i<Traits::nStencilDependencies; ++i)
        stencils.dims[i] = dims[Traits::stencilDependencies[i]];
    stencils.resize(stencils.dims.Product());
    
    typedef std::pair<DiscreteType,OffsetType> IndexOffsetPair;
    std::vector<IndexOffsetPair> offsets;
    offsets.reserve(stencils.size()*HFM::StencilType::nNeigh);
    IndexType updatedIndex;
    auto InsertOffset = [this,pFM,&offsets,&updatedIndex](OffsetType offset, ScalarType w){
        if(w==0.) return;
        IndexType acceptedIndex;
        auto transform = pFM->VisibleOffset(updatedIndex,offset,acceptedIndex);
        if(transform.IsValid()){
            transform.PullVector(offset);
            offsets.push_back({
                stencils.Convert(ShortIndexFromIndex(acceptedIndex)),offset});
        }
    };
    
    for(DiscreteType linearIndex=0; linearIndex<stencils.size(); ++linearIndex){
        updatedIndex = IndexFromShortIndex(stencils.Convert(linearIndex));
        StencilType & stencil = stencils[linearIndex];
        SetStencil(updatedIndex,stencil);
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
    reversedOffsetsSplits.reserve(stencils.size()+1);
    reversedOffsets.reserve(offsets.size());
    
    DiscreteType linearIndex=-1;
    for(const auto & linIndOff : offsets){
        assert(linearIndex<=linIndOff.first);
        for(;linIndOff.first!=linearIndex; ++linearIndex){
            reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());}
        reversedOffsets.push_back(linIndOff.second);
    }
    reversedOffsetsSplits.resize(stencils.size()+1,(DiscreteType)reversedOffsets.size());
}

template<typename Traits> template<typename Dummy> void
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Recomp, Dummy>::
Initialize(const HFM * pFM) {
    shallowStencilQuads.resize(dims.Product());
    
    typedef std::pair<DiscreteType,OffsetType> IndexOffsetPair;
    std::vector<IndexOffsetPair> offsets;
    offsets.reserve(shallowStencilQuads.size()*HFM::StencilType::nNeigh);
    
    IndexType updatedIndex;
    auto InsertOffset = [pFM,&offsets,&updatedIndex](OffsetType offset, ScalarType w){
        if(w==0.) return;
        IndexType acceptedIndex;
        auto transform = pFM->VisibleOffset(updatedIndex,offset,acceptedIndex);
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
        SetStencil(updatedIndex,stencil);
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
    reversedOffsetsSplits.reserve(pFM->values.size()+1);
    reversedOffsets.reserve(offsets.size());
    
    DiscreteType linearIndex=-1;
    for(const auto & linIndOff : offsets){
        assert(linearIndex<=linIndOff.first);
        for(;linIndOff.first!=linearIndex; ++linearIndex){
            reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());}
        reversedOffsets.push_back(linIndOff.second);
    }
    reversedOffsetsSplits.resize(pFM->values.size()+1,(DiscreteType)reversedOffsets.size());
}



#endif /* StencilDataType_h */
