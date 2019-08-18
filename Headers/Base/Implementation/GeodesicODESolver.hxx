// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef GeodesicODESolver_hxx
#define GeodesicODESolver_hxx



template<typename TTraits> GeodesicODESolver<TTraits>::GeodesicODESolver(const HFM & _fm) : GeodesicSolverInterface(_fm) {
    std::vector<IndexType> seedIndices;
    for(const auto & [index,value] : this->fm.seeds) seedIndices.push_back(index);
    targetDistances = LInfDistance(seedIndices, targetTolerance);
}

template<typename Traits> void GeodesicODESolver<Traits>::
Setup(HFMI * that) {
    assert(that!=nullptr);
    IO & io = that->io;
    geodesicStep = io.Get<ScalarType>("geodesicStep",geodesicStep);
    causalityTolerance = io.Get<ScalarType>("geodesicCausalityTolerance",causalityTolerance);
    targetTolerance = (DiscreteType)io.Get<ScalarType>("geodesicTargetTolerance",targetTolerance);
}

template<typename Traits> std::vector<std::vector<typename Traits::PointType> >
GeodesicODESolver<Traits>::Run(HFMI * that, const std::vector<PointType> & tips) {
    assert(that!=nullptr);
    std::vector<std::vector<PointType> > result;
    for(const PointType & tip : tips){
        std::vector<PointType> geodesic;
        geodesic.push_back(tip);
        const bool failed =
            (HFM::DomainType::periodizeUsesBase && !this->fm.dom.PeriodizeNoBase(geodesic.back()).IsValid()) ||
            Run(geodesic);
        if(failed) {Msg() << "Geodesic ODE solver seems to have failed for tip "
            << that->stencil.Param().ReDim( tip ) << ".\n";}
        result.push_back(std::move(geodesic));
    }
    return result;
}


// ----------------- Geodesic flow ----------------


template<typename Traits> struct GeodesicODESolver<Traits>::
FlowAvgType {
    VectorType flow;
    DiscreteType nStationarity;
    DiscreteType targetTolerance;
};

template<typename Traits> auto GeodesicODESolver<Traits>::
GeodesicFlow(const PointType & p, const Array<ShortType,Dimension> & target, FlowCacheType & cache) const -> FlowAvgType {
    const HFM & fm = this->fm;
    IndexType pIndex; // p is in the Voronoi region controlled by this index.
    ScalarType w[Dimension];
    ScalarType sign[Dimension];
    for(int i=0; i<Dimension; ++i) {
        pIndex[i] = (DiscreteType)floor(p[i]);
        w[i] = p[i]-pIndex[i]-0.5;
        sign[i] = w[i]>0 ? 1 : -1;
        w[i] = std::abs(w[i]);
    }
    
    DomainTransformType pIndexTransform;
    if(HFM::DomainType::periodizeUsesBase){
        const auto & dom = fm.dom;
        PointType pIndexCont = dom.PointFromIndex(pIndex);
        pIndexTransform = dom.Periodize(pIndexCont,p);
        pIndex = dom.IndexFromPoint(pIndexCont);
    }
    
    typedef typename FlowCacheType::iterator FlowIteratorType;
    struct WeightedOrientedFlowType {
        FlowIteratorType it;
        ScalarType weight;
        DomainTransformType transform;
    };
    std::array<WeightedOrientedFlowType, (1<<Dimension)> flows;
    DiscreteType
    minTol=std::numeric_limits<DiscreteType>::max(),
    maxStat=0;
    
    // Get the neighbors
    for(int i=0; i< 1<<Dimension; ++i){
        OffsetType offset;
		typedef typename OffsetType::ComponentType OffsetComponentType;
        for(int j=0; j<Dimension; ++j){offset[j] = (i & (1<<j)) ? OffsetComponentType(sign[j]) : 0;}
        IndexType qIndex;
        const auto transform = fm.VisibleOffset(pIndex, offset, qIndex);
        if(!transform.IsValid()) {
            flows[i]={cache.end(),0,transform};
            continue;}
        const DiscreteType qLinearIndex = fm.values.Convert(qIndex);
        minTol=std::min(minTol, (DiscreteType)target[qLinearIndex]);
        
        auto qIt = cache.find(qLinearIndex);
        if(qIt==cache.end()){
            const auto [insertionIt,inserted]=cache.insert({qLinearIndex,{fm.GeodesicFlow(qIndex),0}});
            assert(inserted);
            qIt = insertionIt;
        }
        ScalarType weight = 1;
        for(int j=0; j<Dimension; ++j){weight *= offset[j] ? w[j] : (1-w[j]);}
        
        flows[i] = {qIt,weight,transform};
        maxStat=std::max(maxStat,qIt->second.second);
        ++qIt->second.second;
    }
    
    ScalarType minVal = Traits::Infinity();
    int iMinVal = -1;
    for(int i=0; i< 1<<Dimension; ++i){
        if(flows[i].weight>0){
            const ScalarType val = flows[i].it->second.first.value;
            if(val < minVal){
                minVal=val;
                iMinVal=i;
            }
        }
    }
    
    FlowAvgType result{ VectorType::Constant(0.),maxStat,minTol};
	if(minVal==Traits::Infinity()) return result; // In Wall
    
    const ScalarType minValWidth = flows[iMinVal].it->second.first.width;
//	if(minValWidth==0)  return result; // Not necessarily at seed. Continue
    
    const ScalarType valThreshold =
    minVal+causalityTolerance*minValWidth;
    ScalarType weightSum=0;
    for(const auto & wOFlow : flows){
        const ScalarType weight = wOFlow.weight;
        if(weight<=0) continue;
        const auto & flowData = wOFlow.it->second.first;
        if(flowData.value<=valThreshold){
            weightSum+=weight;
            //            result.value+=weight*flowData.value;
            VectorType inc = weight*flowData.flow;
            wOFlow.transform.PullVector(inc);
            result.flow+=inc;
        }
    }
    if(weightSum>0) {result.flow/=weightSum;} //result.value/=weightSum;
    
    if(HFM::DomainType::periodizeUsesBase){
        pIndexTransform.PullVector(result.flow);}
	
    return result;
}

// ------ Second order Euler scheme -------

template<typename Traits> bool
GeodesicODESolver<Traits>::Run(std::vector<PointType> & geodesic) const {
    if(geodesic.empty()) ExceptionMacro("GeodesicODESolver error : no start point !");
    FlowCacheType flowCache;
    
    const DiscreteType nStationarityMax = DiscreteType(stationarityThreshold/geodesicStep);
    
    std::queue<DiscreteType> targetHistory;
    for(int i=0; i<sqrt(Dimension)/geodesicStep; ++i){
        targetHistory.push(std::numeric_limits<DiscreteType>::max());}
    
    while(true){
        const PointType & p = geodesic.back();
        auto flowData = GeodesicFlow(p,targetDistances,flowCache);
        
        const DiscreteType oldTargetTolerance = targetHistory.front();
        targetHistory.pop();
        targetHistory.push(flowData.targetTolerance);
		if(oldTargetTolerance < flowData.targetTolerance) return false; // Getting away from target.
        
		if(flowData.nStationarity > nStationarityMax) return true; // Stall

		VectorType flow=flowData.flow;
        ScalarType flowNorm = flow.Norm();
		if(flowNorm==0) return false; // At seed
        
        // Do a half step and recompute flow at this point.
        PointType q = p+0.5*geodesicStep*flow/flowNorm;
        
        DomainTransformType qTransform;
        if(HFM::DomainType::periodizeUsesBase){
            [[maybe_unused]] const auto qTransform = this->fm.dom.Periodize(q,p);
            assert(qTransform.IsValid());
        }
        
        flowData = GeodesicFlow(q,targetDistances,flowCache);
        flow=flowData.flow;
        
        if(HFM::DomainType::periodizeUsesBase){
            qTransform.PullVector(flow);}
        
        flowNorm = flow.Norm();
		if(flowNorm==0) return false;
        
        PointType r = p+geodesicStep*flow/flowNorm;
        
        if(HFM::DomainType::periodizeUsesBase){
            [[maybe_unused]] const auto rTransform = this->fm.dom.Periodize(r,p);
            assert(rTransform.IsValid());}
        
        geodesic.push_back(r);
    }
    return true;
}

// ---------- LInf distance ------------
template<typename Traits> auto GeodesicODESolver<Traits>::
LInfDistance(const std::vector<IndexType> & seeds, ShortType upperBound) const -> Array<ShortType,Dimension> {
    const HFM & fm = this->fm;
    // Initialize distance table
    Array<ShortType,Dimension> result;
    result.dims=fm.values.dims;
    result.resize(result.dims.Product(), std::numeric_limits<ShortType>::max());
    std::priority_queue<std::pair<ShortType,IndexType> > queue;
    for(const auto & seed : seeds) {
        result(seed)=0;
        queue.push({0,seed});
    }
    
    // Construct offsets
    std::vector<OffsetType> offsets;
    for(int i=0; i< (1<<2*Dimension); ++i){
        OffsetType offset = OffsetType::Constant(0);
        int j=0;
        for(; j<Dimension; ++j){
            const int k = (i>> (2*j)) & 3;
            if(k==3) break;
            offset[j]=k-1;
        }
        if(j<Dimension) continue;
        if(offset.IsNull()) continue;
        offsets.push_back(offset);
    }
    
    // Run Dijkstra
    while(!queue.empty()){
        const auto & acceptedPair = queue.top();
        const ShortType acceptedValue = -acceptedPair.first; // Minus to get adequate ordering
        const IndexType acceptedIndex =  acceptedPair.second;
        queue.pop();
        if(result(acceptedIndex)!=acceptedValue) continue;
        if(acceptedValue>=upperBound) continue;
        for(const OffsetType & offset : offsets){
            IndexType updatedIndex;
            if(!fm.VisibleOffset(acceptedIndex, offset, updatedIndex).IsValid()) continue;
            ShortType & updatedValue = result(updatedIndex);
            if(updatedValue<=acceptedValue+1) continue;
            updatedValue=acceptedValue+1;
            queue.push({-updatedValue,updatedIndex});
        }
    }
    
    return result;
};


#endif /* GeodesicODESolver_hxx */
