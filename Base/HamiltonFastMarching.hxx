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

template<typename Traits> void
HamiltonFastMarching<Traits>::StencilType::PrintSelf(std::ostream & os) const {
    os << "{";
    if(Traits::nForward) os ExportArrayArrow(forward);
    if(Traits::nSymmetric) os ExportArrayArrow(symmetric);
    if(Traits::nMaxForward) os ExportArrayRecursiveArrow(maxForward, 1);
    if(Traits::nMaxSymmetric) os ExportArrayRecursiveArrow(maxSymmetric, 1);
    os << "}";
}

// --------- Construction -------

template<typename T> HamiltonFastMarching<T>::
HamiltonFastMarching(std::unique_ptr<StencilDataType> _pStencilData):
dom(_pStencilData->dims), pStencilData(std::move(_pStencilData)){
    assert(pStencilData);
    values.dims = pStencilData->dims;
    values.resize(values.dims.ProductOfCoordinates(),Traits::Infinity());
    
    acceptedFlags.dims=values.dims;
    acceptedFlags.resize(values.size(),false);
    
    activeNeighs.dims = values.dims;
    activeNeighs.resize(values.size(),0);
    
    pStencilData->Initialize(this);
};

template<typename T> auto HamiltonFastMarching<T>::
MaxStencilWidth() const -> DiscreteType {
    return std::accumulate(
    pStencilData->reversedOffsets.begin(),pStencilData->reversedOffsets.end(),0,
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
    pStencilData->EraseCache(accepted.linear);
    
    const auto offsets = pStencilData->ReversedOffsets(accepted);
    for(OffsetCRef offset : offsets){
        ConditionalUpdate(accepted.index, offset, top.value);}

    return queue.empty();
}

template<typename T> int HamiltonFastMarching<T>::
PostProcess(IndexCRef acceptedIndex) {
    int result = (sndOrder || !extras.postProcessWithRecompute.empty()) ? Decision::kRecompute : Decision::kAccept;
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
    const auto reversed = VisibleOffset(acceptedIndex,-offset, updated.index);
    if(reversed[Dimension]) return;
    for(int i=0; i<Dimension; ++i)
        if(dom.MayReverse(i) && reversed[i])
            offset[i]*=-1;
    
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

    auto data = pStencilData->UpdateData(updated);
    if(data.quad.minVal==-Traits::Infinity()){
        assert(active==0);
        data.quad.minVal = acceptedValue;
    }

    DiscreteType iNeigh=0;
    for(const auto & diff : data.stencil.forward){
        if(diff.offset==offset){
            assert(!active[iNeigh]);
            active[iNeigh]=1;
            data.quad.Add(acceptedValue,diff.Weight(data.mult));
        }
        ++iNeigh;
    }
    
    for(const auto & diff : data.stencil.symmetric){
        if(diff.offset==offset){
            assert(!active[iNeigh]);
            if(!active[iNeigh+1]){
                active[iNeigh]=1;
                data.quad.Add(acceptedValue,diff.Weight(data.mult));
            }
        }
        if(diff.offset==-offset){
            assert(!active[iNeigh+1]);
            if(!active[iNeigh]){
                active[iNeigh+1]=1;
                data.quad.Add(acceptedValue,diff.Weight(data.mult));
            }
        }
        iNeigh+=2;
    }
    
    for(int iMax=0; iMax<T::nMax; ++iMax){
        for(const auto & diff : data.stencil.maxForward[iMax]){
            if(diff.offset==offset){
                assert(!active[iNeigh]);
                active[iNeigh]=1;
                data.quad.Add(acceptedValue,diff.Weight(data.mult),iMax);
            }
            ++iNeigh;
        }
        
        for(const auto & diff : data.stencil.maxSymmetric[iMax]){
            if(diff.offset==offset){
                assert(!active[iNeigh]);
                if(!active[iNeigh+1]){
                    active[iNeigh]=1;
                    data.quad.Add(acceptedValue,diff.Weight(data.mult), iMax);
                }
            }
            if(diff.offset==-offset){
                assert(!active[iNeigh+1]);
                if(!active[iNeigh]){
                    active[iNeigh+1]=1;
                    data.quad.Add(acceptedValue,diff.Weight(data.mult), iMax);
                }
            }
            iNeigh+=2;
        }
    }
    
    const std::pair<ScalarType,int> result = data.quad.Solve();
    
/*    // Alternatively, only insert in queue if value is strictly decreased.
    ScalarType & val = values[updatedLinearIndex];
    if(result.first>=val) return;*/
       
    values[updated.linear] = result.first;
    
    for(int i=0; i<nMaxBits; ++i) active[nNeigh+i] = result.second & (1<<i);
    queue.push({updated.linear,result.first});
    
}

template<typename T> template<size_t n>
struct HamiltonFastMarching<T>::_QuadType {
    ScalarType minVal=-Traits::Infinity();
    std::pair<ScalarType, int> Solve() const;
    void Add(ScalarType value, ScalarType weight, int);
    void Add(ScalarType value, ScalarType weight);
    PrintSelfMacro(_QuadType);
protected:
    // Note : result is not really needed, except for some assertions.
    struct DataType {ScalarType a=0, b=0, c=-1;
        mutable ScalarType result=Traits::Infinity(); PrintSelfMacro(DataType);};
    std::array<DataType,n> data;
};

template <typename T> template<size_t n> void
HamiltonFastMarching<T>::_QuadType<n>::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(minVal)
    ExportArrayArrow(data)
    << "}";
}

template <typename T> template<size_t n> void
HamiltonFastMarching<T>::_QuadType<n>::DataType::PrintSelf(std::ostream & os) const {
    os << "{" << a << "," << b << "," << c << "," << result << "}";
}

template<typename T> template<size_t n> auto
HamiltonFastMarching<T>::_QuadType<n>::Solve() const -> std::pair<ScalarType, int> {
    for(auto & d : data){
        if(d.a==0) {d.result=Traits::Infinity(); continue;}
        const ScalarType delta2 = d.b*d.b - d.a*d.c;
        const ScalarType delta = sqrt(std::max(0.,delta2));
        d.result = (d.b+delta)/d.a;
//        d.result = std::min(d.result,LargeValue());  // Needed ?
    }
    if(n==0) return {minVal+data[0].result,0};
    int iMin(0);
    ScalarType result = Traits::Infinity();
    for(int i=0; i<data.size(); ++i)
        if(data[i].result<result){
            result=data[i].result;
            iMin=i;
        }
    return {minVal+result,iMin};
}

template<typename T> template<size_t n> void
HamiltonFastMarching<T>::_QuadType<n>::
Add(ScalarType value, ScalarType weight, int i) {
    const ScalarType w=weight, v=value-minVal;
    //  assert(Traits::sndOrder || v>=0);
    const ScalarType wv=w*v, wv2 = (w*v)*v;
    auto & d = data[i];
    //  assert(Traits::sndOrder || d.result>=v);
    d.a+=w;
    d.b+=wv;
    d.c+=wv2;
}

template<typename T> template<size_t n> void
HamiltonFastMarching<T>::_QuadType<n>::
Add(ScalarType value, ScalarType weight) {
    const ScalarType w=weight, v=value-minVal;
    //  assert(Traits::sndOrder || v>=0);
    const ScalarType wv=w*v, wv2 = (w*v)*v;
    for(auto & d : data){
//        assert(Traits::sndOrder || d.result>=v);
        d.a+=w;
        d.b+=wv;
        d.c+=wv2;
    }
}



// ----------------- Boundary conditions -------------------

template<typename T> auto HamiltonFastMarching<T>::
VisibleOffset(const IndexType & acceptedIndex, const OffsetType & offset, IndexType & updatedIndex) const -> ReverseFlag {
    for(int i=0; i<Dimension; ++i) {updatedIndex[i]=acceptedIndex[i]+offset[i];}
    std::bitset<Dimension+1> result = dom.Periodize(updatedIndex);
    for(ExtraAlgorithmInterface * p : extras.visible){
        result[Dimension] = result[Dimension] || !p->Visible(acceptedIndex,offset,updatedIndex);}
    return result;
}

// ---------- Recompute ----------

template<typename T> auto HamiltonFastMarching<T>::
Recompute(const IndexType & updatedIndex, DiscreteFlowType & discreteFlow) const -> RecomputeType {
    for(ExtraAlgorithmInterface * p : extras.beforeRecompute) p->BeforeRecompute(updatedIndex);
    assert(discreteFlow.empty());
    const DiscreteType updatedLinearIndex = values.Convert(updatedIndex);
    const ActiveNeighFlagType active = activeNeighs[updatedLinearIndex];
//    if(active.none()) return values[updatedLinearIndex]; // Handled below
    
    const auto & data = pStencilData->RecomputeData(updatedIndex);
    
    std::bitset<nActiveNeigh> sndOrderNeighs;
    auto PushValueDiff = [this,&updatedIndex, &discreteFlow, &sndOrderNeighs](OffsetType offset){
        IndexType acceptedIndex;
        for(int i=0; i<Dimension; ++i) {acceptedIndex[i] = updatedIndex[i]+offset[i];}
        const auto & flipped = dom.Periodize(acceptedIndex);
        assert(!flipped[Dimension]);
        
        //Note that, in the beginning, we store values, and not weights, in discreteFlow.
        const ScalarType acceptedValue = values(acceptedIndex);
        discreteFlow.push_back({offset,acceptedValue});
        
        if(!sndOrder) return;
        OffsetType offset2;
        for(int i=0; i<Dimension; ++i){offset2[i] = flipped[i] ? -offset[i] : offset[i];}
        IndexType acceptedIndex2;
        if(VisibleOffset(acceptedIndex, offset2, acceptedIndex2)[Dimension]) return;
        
        const DiscreteType acceptedLinearIndex2 = values.Convert(acceptedIndex2);
        if(!acceptedFlags[acceptedLinearIndex2]) return;
        const ScalarType acceptedValue2 = values(acceptedIndex2);
        if(acceptedValue2>acceptedValue) return;
        discreteFlow.back().weight = (4./3.)*acceptedValue - (1./3.)*acceptedValue2;
        sndOrderNeighs[discreteFlow.size()-1]=1;
    };
    
    CappedVector<ScalarType, nActiveNeigh> weights;
    int iNeigh=0;
    for(const auto & diff : data.stencil.forward){
        if(active[iNeigh]){
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff(diff.offset);
        }
        ++iNeigh;
    }
    for(const auto & diff : data.stencil.symmetric){
        if(active[iNeigh]){
            assert(!active[iNeigh+1]);
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff( diff.offset);
        }
        if(active[iNeigh+1]){
            assert(!active[iNeigh]);
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff(-diff.offset);
        }
        iNeigh+=2;
    }
    
    int iMax=0;
    for(int i=0; i<nMaxBits; ++i){
        if(active[nNeigh+i]) iMax |= (1<<i);}
    iNeigh+=iMax*(T::nMaxForward+2*T::nMaxSymmetric);
    const auto & maxForward = data.stencil.maxForward[iMax];
    const auto & maxSymmetric = data.stencil.maxSymmetric[iMax];
    
    for(const auto & diff : maxForward){
        if(active[iNeigh]){
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff(diff.offset);
        }
        ++iNeigh;
    }
    for(const auto & diff : maxSymmetric){
        if(active[iNeigh]){
            assert(!active[iNeigh+1]);
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff(diff.offset);
        }
        if(active[iNeigh+1]){
            assert(!active[iNeigh]);
            weights.push_back(diff.Weight(data.mult));
            PushValueDiff(-diff.offset);
        }
        iNeigh+=2;
    }
    
    _QuadType<1> quad; quad.minVal = Traits::Infinity();
    for(const auto & offsetValue : discreteFlow){
        quad.minVal=std::min(quad.minVal, offsetValue.weight);}
    
    for(int i=0; i<discreteFlow.size(); ++i){
        quad.Add(discreteFlow[i].weight,
                 (sndOrder && sndOrderNeighs[i]) ?
                 (9./4.)* weights[i] : weights[i]);
    }

/*    if(discreteFlow.empty()){
        std::cout
        ExportVarArrow(iNeigh)
        ExportArrayArrow(weights)
        ExportVarArrow(active)
        ExportVarArrow(iMax)
        << "\n";
    }*/

    if(discreteFlow.empty()) return {values[updatedLinearIndex],0.};
    RecomputeType result;
    result.value = quad.Solve().first;
    result.width = 0.;
    
    for(int i=0; i<discreteFlow.size(); ++i){
        // Difference should already be non-negative, by construction, without sndOrder.
        const ScalarType weightPosDiff = weights[i]*std::max(0.,result.value-discreteFlow[i].weight);
        discreteFlow[i].weight = weightPosDiff;
        result.width+=weightPosDiff;        
        if(sndOrder && sndOrderNeighs[i]) discreteFlow[i].weight *= 3./2.;
    }
    const ScalarType weightsSum = std::accumulate(weights.begin(), weights.end(), 0.);
    if(weightsSum>0) {result.width/=weightsSum;}
    return result;
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
