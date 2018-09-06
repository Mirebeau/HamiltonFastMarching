//
//  DynamicFactoring.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/08/2018.
//

#ifndef DynamicFactoring_hxx
#define DynamicFactoring_hxx

template<> char const* enumStrings<FactoringPointChoice>::data[] = {"Key", "Current", "Both"};

template<typename T> void
DynamicFactoring<T>::ElementaryGuess::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(guess)
    ExportVarArrow(base)
    ExportVarArrow(weight)
    << "}";
}

template<typename T> auto
DynamicFactoring<T>::
Correction(const OffsetType & offset, bool snd) const -> ScalarType {
    /*std::cout
    ExportArrayArrow(guesses)
    ExportVarArrow(offset)
    << "\n";*/
//    if(true || pointChoice!=FactoringPointChoice::Both){
        return std::accumulate(guesses.begin(),guesses.end(),0.,
                               [&](ScalarType a, const ElementaryGuess & b)->ScalarType{
                                   return a + b.Correction(offset,snd);});
    /*} else {
        const int n2 = guesses.size();
        assert(n2%2==0);
        const int n = n2/2;
        return
        std::accumulate(guesses.begin(),guesses.begin()+n,0.,
                        [&](ScalarType a, const ElementaryGuess & b)->ScalarType{
                            return a + b.Correction(offset,snd,2.);}) +
        std::accumulate(guesses.begin()+n,guesses.end(),0.,
                        [&](ScalarType a, const ElementaryGuess & b)->ScalarType{
                            return a + b.Correction(offset,snd,0.);});
        
    }*/
}

template<typename T> bool // Returns wether dynamic factoring is active at this point
DynamicFactoring<T>::
SetIndex(IndexCRef index, const DiscreteFlowType & flow) {
    assert(!factoringRegion.empty());
    const DiscreteType linearIndex = factoringRegion.Convert(index);
    if(!factoringRegion[linearIndex]) return false;

//    std::cout ExportVarArrow(index) << "\n";
    const FullIndexCRef full{index,linearIndex};
    MakeFactor(full, flow);
    return MakeGuess(full);
}

template<typename T> void
DynamicFactoring<T>::
MakeFactor(FullIndexCRef updated, const DiscreteFlowType & flow){
    if(factoringDone[updated.linear]) return; // Work already done
    factoringDone[updated.linear]=true;
    
    currentNeighbors.clear();
    currentNeighbors2.clear();
    
    
    const ScalarType wSum = std::accumulate(flow.begin(), flow.end(), 0.,
         [](ScalarType a, const DiscreteFlowElement & b)->ScalarType{return a+b.weight;});
    assert(wSum>0);
    
    for(const auto & offsetWeight : flow){
        const OffsetType & offset = offsetWeight.offset;
        const ScalarType & weight = offsetWeight.weight/wSum;
        
        IndexType neigh = updated.index+IndexDiff::CastCoordinates(offset);
        assert(pFM!=nullptr);
        const auto transform = pFM->dom.Periodize(neigh,updated.index);
        assert(transform.IsValid());
        const DiscreteType linearNeigh = factoringRegion.Convert(neigh);
        
        if(factoringKeypoints.find(linearNeigh)!=factoringKeypoints.end()){
            currentNeighbors.push_back({offset,weight});
        } else {
            const auto eq = factoringNeighbors.equal_range(linearNeigh);
            for(auto it = eq.first; it!=eq.second; ++it){
                const auto & offsetWeight = it->second;
                currentNeighbors.push_back({offsetWeight.first + offset, offsetWeight.second * weight});
                // TODO : remove some of these, based on offset radius ?
            }
        }
    }
    
#ifdef Debug
    const ScalarType weightSum = std::accumulate(currentNeighbors.begin(),currentNeighbors.end(),0.,[](ScalarType a, const std::pair<OffsetType,ScalarType> & b) -> ScalarType {return a+b.second;});
    assert(weightSum<1.001);
#endif
    
    // Merge duplicates and insert
    std::sort(currentNeighbors.begin(),currentNeighbors.end());
    for(const auto & offsetWeight : currentNeighbors){
        if(currentNeighbors2.empty() || currentNeighbors2.back().first != offsetWeight.first){
            currentNeighbors2.push_back(offsetWeight);
        } else {
            currentNeighbors2.back().second+=offsetWeight.second;
        }
    }
    
    const auto it = factoringNeighbors.upper_bound(updated.linear);
    for(const auto & data : currentNeighbors2){
        factoringNeighbors.insert(it,{updated.linear,data});}
}

template<typename T> bool
DynamicFactoring<T>::
MakeGuess(FullIndexCRef updated) {
    assert(factoringRegion[updated.linear]);
    assert(factoringDone[updated.linear]);
    guesses.clear();
    const auto rg = factoringNeighbors.equal_range(updated.linear);
    if(rg.first==rg.second) return false;
    
    if(pointChoice==FactoringPointChoice::Current || pointChoice==FactoringPointChoice::Both){
        const DistanceGuess guess = pFM->stencilData.GetGuess(updated.index);
        for(auto it = rg.first; it!=rg.second; ++it){
            const auto & offsetWeight = it->second;
            const DomainTransformType transform{}; // Identity transform
            guesses.push_back({guess,transform,offsetWeight.first,
                offsetWeight.second *(pointChoice==FactoringPointChoice::Both ? 0.5 : 1.)});
        }
    }
    if(pointChoice==FactoringPointChoice::Key || pointChoice==FactoringPointChoice::Both){
        for(auto it = rg.first; it!=rg.second; ++it){
            const auto & offsetWeight = it->second;
            IndexType neigh = updated.index+IndexDiff::CastCoordinates(offsetWeight.first);
            const DomainTransformType transform = pFM->dom.Periodize(neigh,updated.index);
            const auto keyIt = factoringKeypoints.find(factoringDone.Convert(neigh));
            assert(keyIt!=factoringKeypoints.end());
            guesses.push_back({keyIt->second,transform,offsetWeight.first,
                offsetWeight.second *(pointChoice==FactoringPointChoice::Both ? 0.5 : 1.)});
        }
    }
    return true;
}

template<typename T> bool
DynamicFactoring<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(!io.template Get<ScalarType>("useFactoring",0.,1)) return false;
        
    factoringRadius = io.template Get<ScalarType>("factoringRadius",factoringRadius);
    pointChoice = enumFromString<FactoringPointChoice>(io.GetString("factoringPointChoice",
                                                                    enumToRealString(pointChoice)));
    if(0>(int)pointChoice) {ExceptionMacro("Dynamic factoring error: unrecognized factoringPointChoice.");}
    
    pFM = that->pFM.get();
    const auto & stencilData = pFM->stencilData;
    const auto & param = pFM->stencilData.Param();
    const auto & dom = pFM->dom;
    
    // Setup arrays
    const auto & dims = stencilData.dims;
    const size_t size = dims.ProductOfCoordinates();
    factoringRegion.dims = dims; factoringRegion.resize(size,false);
    factoringDone.dims = dims; factoringDone.resize(size,false);
    const auto & arr = factoringDone;
    
    // ------ Setup some edges to be used in Dijkstra methods at initialization ------
    { // Setup the edges
        Array<int, Dimension> dummy; dummy.dims = IndexType::Constant(3);
        int nEdges = dummy.dims.ProductOfCoordinates();
        for(int i=0; i<nEdges; ++i){
            const IndexDiff offset = dummy.Convert(i)-IndexType::Constant(1);
            const ScalarType norm = VectorType::CastCoordinates(offset).Norm();
            if(norm>0){edges.push_back({offset,norm});}
        }
    }
    
    // Setup keypoints
    const std::vector<DiscreteType> keypoints = SelectKeypoints(that);
    if(io.template Get<ScalarType>("exportFactoringKeypoints",0.,2)){
        std::vector<PointType> pts;
        for(DiscreteType i : keypoints) {
            pts.push_back(param.ReDim(dom.PointFromIndex(arr.Convert(i))));}
        io.SetVector("factoringKeypoints",pts);
    }
    
    std::priority_queue<std::pair<ScalarType,DiscreteType> > queue;
    for(DiscreteType i : keypoints){
        queue.push({0.,i});
        factoringKeypoints.insert({i,stencilData.GetGuess(arr.Convert(i))});
    }
    
    // ------ Setup dynamic factoring region, by running a Dijkstra, approximating euclidean distance ------
    while(!queue.empty()){
        const auto top = queue.top();
        const ScalarType value = - top.first;
        const DiscreteType current = top.second;
        queue.pop();
        if(factoringRegion[current]) continue;
        factoringRegion[current]=true;
        for(const auto & e : edges){
            const ScalarType neighVal = value+e.second;
            if(neighVal>factoringRadius) continue;
            const IndexType index = arr.Convert(current);
            IndexType neigh = index+e.first;
            if(dom.Periodize(neigh,index).IsValid()){
                queue.push({-neighVal,arr.Convert(neigh)});}
        }
    }
    
    if(io.template Get<ScalarType>("exportFactoringRegion",0.,2.)){
        io.SetArray("factoringRegion",factoringRegion.template Cast<ScalarType>());}
    return true;
}

template<typename T> auto
DynamicFactoring<T>::
SelectKeypoints(HFMI * that) -> std::vector<DiscreteType> {
    std::vector<DiscreteType> result;
    
    IO & io = that->io;
    const auto & param = pFM->stencilData.Param();
    const auto & dom = pFM->dom;

    Array<bool, Dimension> & done = factoringDone;
    assert(!done.empty());
    
    auto LinearIndexFromPhysicalPoint = [&](const PointType & p){
        IndexType index = dom.IndexFromPoint(param.ADim(p));
        if(!dom.PeriodizeNoBase(index).IsValid()){
            ExceptionMacro("Dynamic factoring error : input point " << p << " is out of range");}
        return done.Convert(index);
    };
    
    if(io.HasField("factoringKeypoints")){
        const std::vector<PointType> pts = io.template GetVector<PointType>("factoringKeypoints");
        for(const PointType & p : pts){result.push_back(LinearIndexFromPhysicalPoint(p));}
        return result;
    }
    
    const std::vector<PointType> seeds = io.template GetVector<PointType>("seeds");
    for(const PointType & p : seeds){result.push_back(LinearIndexFromPhysicalPoint(p));}

    // ------ Extract keypoints from wall data -----
    
    if(!io.HasField("walls")) return result;
    const auto pWalls = that->template GetIntegralField<bool>("walls");
    pWalls->CheckDims(done.dims);
    assert(!edges.empty());
    
    const ScalarType factoringWallExtractionRadius = io.template Get<ScalarType>("factoringWallExtractionRadius",factoringRadius,2);
    
    // Get the wall boundary, which are the candidate keypoints
    std::set<DiscreteType> wallBoundary;
    for(DiscreteType i=0; i<done.size(); ++i){
        const IndexType index = done.Convert(i);
        if(!(*pWalls)(index)) continue;
        for(const auto & e : edges){
            IndexType neigh = index +e.first;
            if(!dom.Periodize(neigh,index).IsValid()) continue;
            if((*pWalls)(neigh)) continue;
            wallBoundary.insert(done.Convert(neigh));
        }
    }
    std::priority_queue<std::tuple<ScalarType,DiscreteType,DiscreteType> > queue; //-value,current,source
    for(DiscreteType i : wallBoundary){ queue.push({0.,i,i});}
    
    // Find the singular points of the wall boundary, by growing balls around them
    std::map<DiscreteType,DiscreteType> score;
    while(!queue.empty()){
        const auto top = queue.top();
        const ScalarType value = - std::get<0>(top);
        const DiscreteType current = std::get<1>(top), source = std::get<2>(top);
        queue.pop();
        if(done[current]) continue;
        done[current]=true;
        if(source!=current){
            //            std::cout ExportVarArrow(done.Convert(source))  ExportVarArrow(done.Convert(current)) << "\n";
            const auto it=score.find(source);
            if(it!=score.end()) {++(it->second);}
            else {score.insert({source,1});}
        }
        for(const auto & e : edges){
            const ScalarType neighVal = value+e.second;
            if(neighVal>factoringWallExtractionRadius) continue;
            const IndexType index = done.Convert(current);
            IndexType neigh = index+e.first;
            if(dom.Periodize(neigh,index).IsValid() && !(*pWalls)(neigh)){
                queue.push({-neighVal,done.Convert(neigh),source});}
        }
    }
    
    
    const ScalarType factoringWallExtractionThreshold = io.template Get<ScalarType>("factoringWallExtractionThreshold", factoringRadius*sqrt(factoringRadius), 2);
    for(const auto & indexScore : score){
        if(indexScore.second>factoringWallExtractionThreshold && !OnBoundary(done.Convert(indexScore.first), dom)){
            result.push_back(indexScore.first); } }
    
    // Eliminate erroneous keypoints, lying on the boundary
    // Alternatively, regard the boundary as an obstacle.
    
/*    std::cout
    ExportVarArrow(factoringWallExtractionRadius)
    ExportVarArrow(factoringWallExtractionThreshold)
    ExportArrayArrow(score) << "\n";*/
    
    std::fill(done.begin(),done.end(),false); // reset this array for future use
    return result;
}

template<typename T> bool
DynamicFactoring<T>::
OnBoundary(IndexCRef index, const DomainType & dom) const {
    for(int i=0; i<Dimension; ++i){
        for(int eps=-1; eps<=1; eps+=2){
            IndexType neigh=index;
            neigh[i]+=eps;
            if(!dom.Periodize(neigh,index).IsValid()) return true;
        }
    }
    return false;
}

#endif /* DynamicFactoring_hxx */
