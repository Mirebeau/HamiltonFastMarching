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

// --- Update data ---
/*
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<true, Dummy>::
UpdateData(FullIndexCRef full) -> UpdateDataType {
    const bool found = shallowMultQuads.find(full.linear);
    auto & multQuad = shallowMultQuads[full.linear];
    if(!found) { multQuad.first = (*pMultSource)(full.index);}
    const StencilType & stencil = stencils(ShortIndexFromIndex(full.index));
    
    return UpdateDataType{stencil,multQuad.first,multQuad.second};
}*/

/*
template<typename Traits> template<typename Dummy> auto
HamiltonFastMarching<Traits>::_StencilDataType<false, Dummy>::
UpdateData(FullIndexCRef full) -> UpdateDataType {
    const bool found = shallowStencilQuads.find(full.linear);
    auto & stencilQuad = shallowStencilQuads[full.linear];
    if(!found) {SetStencil(full.index,stencilQuad.first);}
    
    return UpdateDataType{stencilQuad.first,MultiplierType{},stencilQuad.second};
}
*/
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
/*    std::cout
    ExportVarArrow(index)
    ExportVarArrow(i)
    ExportArrayArrow(reversedOffsets)
    ExportArrayArrow(reversedOffsetsSplits)
    << "\n";*/
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
        /*
        for(const DifferenceType & diff : stencil.forward)
            InsertOffset( diff.offset, diff.baseWeight);
        for(const DifferenceType & diff : stencil.symmetric){
            InsertOffset( diff.offset, diff.baseWeight);
            InsertOffset(-diff.offset, diff.baseWeight);
        }
        for(const auto & diffs : stencil.maxForward)
            for(const auto & diff : diffs)
                InsertOffset( diff.offset, diff.baseWeight);
        for(const auto & diffs : stencil.maxSymmetric)
            for(const auto & diff : diffs){
                InsertOffset( diff.offset, diff.baseWeight);
                InsertOffset(-diff.offset, diff.baseWeight);
            }
         */
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
    
/*    std::cout
    ExportArrayArrow(reversedOffsets)
    ExportArrayArrow(reversedOffsetsSplits)
    << "\n";*/
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
/*        for(const DifferenceType & diff : stencil.forward)
            InsertOffset( diff.offset, diff.baseWeight);
        for(const DifferenceType & diff : stencil.symmetric){
            InsertOffset( diff.offset, diff.baseWeight);
            InsertOffset(-diff.offset, diff.baseWeight);
        }*/
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
    //    [](const IndexOffsetPair & a, const IndexOffsetPair & b){return a.first<b.first;}
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
    
    directOffsetsSplits.reserve(size);
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
    for(const auto indexOffset : targets){
        reversedOffsets.push_back(indexOffset.second);
        if(last!=indexOffset.first){
            reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());
            last = indexOffset.first;
        }
    }
    reversedOffsetsSplits.push_back((DiscreteType)reversedOffsets.size());
}

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
    
    const auto valAct = HopfLaxUpdate(updatedIndex,offsetVal);
    const ScalarType oldValue = pFM->values[full.linear];
    
    if(valAct.first>=oldValue) return oldValue;
    
    active.sectorIndex = act[valAct.second];
    return valAct.first;
}

template<typename Traits> template<typename Dummy> template<typename F> auto
HamiltonFastMarching<Traits>::_StencilDataType<SSP::Lag2, Dummy>::
HopfLaxRecompute(const F & f, IndexCRef index, ActiveNeighFlagType active, DiscreteFlowType & discreteFlow) -> RecomputeType {
    assert(!active.none());
    StencilType stencil;
    SetStencil(index, stencil);
    
    // Get data of neighbor 0
    int snd0 = 1;
    const OffsetType offset0 = stencil.Sector(active.sectorIndex,0);
    ScalarType val0 = f(offset0,snd0);
    
    // Get data of neighbor 1
    int snd1 = snd0==0 ? 0 : 1;
    const OffsetType offset1 = stencil.Sector(active.sectorIndex,1);
    const ScalarType val1 = f(offset1,snd1);

    if(snd1==0 && snd0==1){
        snd0=0;
        val0 = f(offset0,snd0);
        assert(snd0==0);
    }
    
    assert(snd0!=-1 || snd1!=-1);
    int snd = snd0==-1 ? snd1 : snd0; // snd order used ?
    assert(snd!=-1);
    assert((snd0==-1 || snd0==snd) && (snd1==-1 || snd1==snd));
    
    assert(discreteFlow.empty());
    if(snd0!=-1) {discreteFlow.push_back({offset0,snd==0 ? val0 : (1.5*val0)});}
    if(snd1!=-1) {discreteFlow.push_back({offset1,snd==0 ? val1 : (1.5*val1)});}
    
    RecomputeType result = HopfLaxRecompute(index,discreteFlow);
    if(snd){result.width/=1.5; result.value/=1.5;}
    
    return result;
};


#endif /* StencilDataType_h */
