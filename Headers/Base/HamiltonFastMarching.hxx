// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HamiltonFastMarching_hxx
#define HamiltonFastMarching_hxx

// ----- Printing some structures ------

template<typename Traits> void
HamiltonFastMarching<Traits>::FlowDataType::PrintSelf(std::ostream & os) const {
    os << "{" << flow << "," << value << "," << width << "}";
}
/*
template<typename Traits> void
HamiltonFastMarching<Traits>::StencilType::PrintSelf(std::ostream & os) const {
    os << "{";
    if(Traits::nForward) os ExportArrayArrow(forward);
    if(Traits::nSymmetric) os ExportArrayArrow(symmetric);
    if(Traits::nMaxForward) os ExportArrayRecursiveArrow(maxForward, 1);
    if(Traits::nMaxSymmetric) os ExportArrayRecursiveArrow(maxSymmetric, 1);
    os << "}";
}*/

// --------- Construction -------

template<typename T> HamiltonFastMarching<T>::
HamiltonFastMarching(StencilDataType & _stencilData):
dom(_stencilData.dims), stencilData(_stencilData){
    values.dims = stencilData.dims;
    values.resize(values.dims.ProductOfCoordinates(),Traits::Infinity());
    
    acceptedFlags.dims=values.dims;
    acceptedFlags.resize(values.size(),false);
    
    activeNeighs.dims = values.dims;
    activeNeighs.resize(values.size());
    
    stencilData.Initialize(this);
};

template<typename T> auto HamiltonFastMarching<T>::
MaxStencilWidth() const -> DiscreteType {
    return std::accumulate(
    stencilData.reversedOffsets.begin(),stencilData.reversedOffsets.end(),0,
    [](DiscreteType a, const OffsetType & offset)->DiscreteType {
    return std::max(a,std::accumulate(offset.begin(),offset.end(),0,
    [](DiscreteType b, DiscreteType c)->DiscreteType {return b+std::abs(c);}
                                      ));});
}


// ------------ Running Fast Marching -------------


template<typename T>
struct HamiltonFastMarching<T>::QueueElement {
    DiscreteType linearIndex;
    ScalarType value;
    bool operator < (const QueueElement & other) const {return value>other.value;}
};


template<typename T>
void HamiltonFastMarching<T>::Run(){
    RunInit();
    do {} while(!RunOnce());
}

template<typename T>
void HamiltonFastMarching<T>::RunInit(){
    assert(values.CheckDims());
    assert(values.size()>0);
    
    if(acceptedFlags.empty()){
        acceptedFlags.dims = values.dims;
        acceptedFlags.resize(values.size(),false);
    } else {
        assert(acceptedFlags.dims==values.dims);
        assert(acceptedFlags.CheckDims());
    }
    
    for(const auto & seed : seeds){
        const DiscreteType linearIndex = values.Convert(seed.first);
        values[linearIndex] = seed.second;
//        acceptedFlags[linearIndex] = true;
        queue.push({linearIndex,seed.second});
    }
}

template<typename T>
bool HamiltonFastMarching<T>::RunOnce(){
    const QueueElement top = queue.top();
    queue.pop();
    
    if(acceptedFlags[top.linearIndex]) return queue.empty();
    const FullIndexType accepted = {values.Convert(top.linearIndex),top.linearIndex};
    int dec = PostProcess(accepted.index);
    if(dec & Decision::kRecompute){
        DiscreteFlowType flow;
        const RecomputeType rec = Recompute(accepted.index, flow);
        dec|=PostProcessWithRecompute(accepted.index, rec, flow);
    }
    if(dec & Decision::kTerminate) return true;
    if(dec & Decision::kContinue) return queue.empty();
    
    acceptedFlags[accepted.linear]=true;
    stencilData.EraseCache(accepted.linear);
    
    const auto offsets = stencilData.ReversedOffsets(accepted);
    for(OffsetCRef offset : offsets){
        ConditionalUpdate(accepted.index, offset, top.value);}

    return queue.empty();
}

template<typename T> int HamiltonFastMarching<T>::
PostProcess(IndexCRef acceptedIndex) {
    int result = (sndOrder || dynamicFactoring!=nullptr || !extras.postProcessWithRecompute.empty()) ?
    Decision::kRecompute : Decision::kAccept;
    for(ExtraAlgorithmInterface * p : extras.postProcess) result|=p->PostProcess(acceptedIndex);
    return result;
}

template<typename T> int HamiltonFastMarching<T>::
PostProcessWithRecompute(IndexCRef acceptedIndex, const RecomputeType & rec, const DiscreteFlowType & flow){
    values(acceptedIndex)=rec.value;
    int result = Decision::kAccept;
    for(ExtraAlgorithmInterface * p : extras.postProcessWithRecompute)
        result|=p->PostProcessWithRecompute(acceptedIndex, rec, flow);
    return result;
}

template<typename T> void
HamiltonFastMarching<T>::ConditionalUpdate(IndexCRef acceptedIndex,
                                           OffsetType offset,
                                           ScalarType acceptedValue){
    FullIndexType updated;
    const auto transform = VisibleOffset(acceptedIndex,-offset, updated.index);
    if(!transform.IsValid()) return;
    transform.PullVector(offset);    
    updated.linear = values.Convert(updated.index);
    if(acceptedFlags[updated.linear]) return;

    // Next line forbids updating of seeds or given data.
    // Can be specialized to allow e.g. for sequential computation of Voronoi diagrams.
    if(activeNeighs[updated.linear].none() && values[updated.linear]!=Traits::Infinity()) return;
    
    if(values[updated.linear]<=acceptedValue) return;
    Update(updated, offset, acceptedValue); // also pushes in queue
}


template<typename T> void HamiltonFastMarching<T>::
Update(FullIndexCRef updated, OffsetCRef offset, ScalarType acceptedValue){
    auto & active = activeNeighs[updated.linear];

    
/*    // Alternatively, only insert in queue if value is strictly decreased.
    ScalarType & val = values[updatedLinearIndex];
    if(result.first>=val) return;*/
    
    const ScalarType updatedValue =
    stencilData.HopfLaxUpdate(updated,offset,acceptedValue,active);
    values[updated.linear] = updatedValue;
    queue.push({updated.linear,updatedValue});
    
}

// ----------------- Boundary conditions -------------------

template<typename T> auto HamiltonFastMarching<T>::
VisibleOffset(const IndexType & acceptedIndex, const OffsetType & offset, IndexType & updatedIndex) const -> DomainTransformType {
    updatedIndex=acceptedIndex+IndexDiff::CastCoordinates(offset);
    DomainTransformType result = dom.Periodize(updatedIndex,acceptedIndex);
    if(!result.IsValid()) return result;
    for(ExtraAlgorithmInterface * p : extras.visible){
        if(!p->Visible(acceptedIndex,offset,updatedIndex))
            result.Invalidate();
            }
    return result;
}

// ---------- Recompute ----------

template<typename T> auto HamiltonFastMarching<T>::
Recompute(IndexCRef updatedIndex, DiscreteFlowType & discreteFlow) const -> RecomputeType {
    assert(discreteFlow.empty());
    for(ExtraAlgorithmInterface * p : extras.beforeRecompute) {p->BeforeRecompute(updatedIndex);}
    const DiscreteType updatedLinearIndex = values.Convert(updatedIndex);
    const ActiveNeighFlagType active = activeNeighs[updatedLinearIndex];
    if(active.none()) return {values[updatedLinearIndex],0.}; // Handled below
    
    auto GetValue = [this,&updatedIndex](OffsetType offset, int & snd) -> ScalarType {
        //snd code : -1 -> invalid, 0 -> first order, 1 -> sndOrder.
        
        IndexType acceptedIndex = updatedIndex+IndexDiff::CastCoordinates(offset);
        const auto transform = dom.Periodize(acceptedIndex,updatedIndex);
        if(!transform.IsValid()) {snd=-1; return -Traits::Infinity();}
        const ScalarType acceptedValue = values(acceptedIndex);
        
        const bool trySnd = sndOrder && snd;
        while(true){
            if(!trySnd) break;
            OffsetType offset2 = offset;
            transform.PullVector(offset2);
            IndexType acceptedIndex2;
            if(!VisibleOffset(acceptedIndex, offset2, acceptedIndex2).IsValid()) break;
            
            const DiscreteType acceptedLinearIndex2 = values.Convert(acceptedIndex2);
            if(!acceptedFlags[acceptedLinearIndex2]) break;
            const ScalarType acceptedValue2 = values(acceptedIndex2);
            if(acceptedValue2>acceptedValue) break;
            snd=1;
            return (4./3.)*acceptedValue - (1./3.)*acceptedValue2;
        }
        snd=0;
        return acceptedValue;
    };
    
    const RecomputeType & rec = stencilData.HopfLaxRecompute(GetValue,updatedIndex,active,discreteFlow);
    
    if(dynamicFactoring==nullptr || !dynamicFactoring->SetIndex(updatedIndex,discreteFlow)){
        return rec;}
    
    // TODO : Choose wether sndOrder and factoring should be enable together. (Seems not.)
    
    auto GetValueCorr = [this,&updatedIndex](OffsetType offset, int & snd) -> ScalarType {
        //snd code : -1 -> invalid, 0 -> first order, 1 -> sndOrder.
        
        IndexType acceptedIndex = updatedIndex+IndexDiff::CastCoordinates(offset);
        const auto transform = dom.Periodize(acceptedIndex,updatedIndex);
        if(!transform.IsValid()) {snd=-1; return -Traits::Infinity();}
        const ScalarType acceptedValue = values(acceptedIndex);
        
        const bool trySnd = sndOrder && snd;
        while(false && true){
            if(!trySnd) break;
            OffsetType offset2 = offset;
            transform.PullVector(offset2);
            IndexType acceptedIndex2;
            if(!VisibleOffset(acceptedIndex, offset2, acceptedIndex2).IsValid()) break;
            const DiscreteType acceptedLinearIndex2 = values.Convert(acceptedIndex2);
            if(!acceptedFlags[acceptedLinearIndex2]) break;
            const ScalarType acceptedValue2 = values(acceptedIndex2);
            if(acceptedValue2>acceptedValue) break;
            snd=1;
            return (4./3.)*acceptedValue - (1./3.)*acceptedValue2
            + dynamicFactoring->Correction(offset,true);
        }
        snd=0;
        return acceptedValue + dynamicFactoring->Correction(offset,false);
    };
    
    
/*    if(updatedIndex==IndexType{14,12}){
        discreteFlow.clear();
        const RecomputeType rec2 = stencilData.HopfLaxRecompute(GetValueCorr,updatedIndex,active,discreteFlow);
        int snd0 = 0;
        std::cout << "Recomputing "
        ExportVarArrow(updatedIndex)
        ExportVarArrow(GetValue(OffsetType{-1,0},snd0))
//        ExportVarArrow(GetValue(OffsetType{0,-1},snd0))
        ExportVarArrow(GetValueCorr(OffsetType{-1,0},snd0))
//        ExportVarArrow(GetValueCorr(OffsetType{0,-1},snd0))
        ExportVarArrow(rec.value)
        ExportVarArrow(rec2.value)
        << "\n";
        
    }*/
    
    discreteFlow.clear();
    return stencilData.HopfLaxRecompute(GetValueCorr,updatedIndex,active,discreteFlow);
    
}

template<typename Traits> auto HamiltonFastMarching<Traits>::
GeodesicFlow(const IndexType & index) const -> FlowDataType {
    DiscreteFlowType discreteFlow;
    const RecomputeType & rec = Recompute(index,discreteFlow);
    FlowDataType result;
    result.value =  rec.value;
    result.width = rec.width;
    
    result.flow.fill(0.);
    for(const auto & offsetWeight : discreteFlow){
        const OffsetType & offset = offsetWeight.offset;
        const ScalarType & weight = offsetWeight.weight;
        assert(weight>=0);
        for(int i=0; i<Dimension; ++i){
            result.flow[i]+=weight*ScalarType(offset[i]);}
    }
    return result;
}


#endif /* HamiltonFastMarching_hxx */
