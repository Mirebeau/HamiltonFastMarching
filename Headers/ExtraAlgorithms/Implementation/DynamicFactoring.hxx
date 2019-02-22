//
//  DynamicFactoring.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/08/2018.
//

#ifndef Factoring_hxx
#define Factoring_hxx

template<> char const* enumStrings<FactoringPointChoice>::data[] = {"Key", "Current", "Both"};
template<> char const* enumStrings<FactoringMethod>::data[] = {"None", "Static", "Dynamic"};

template<typename T> void
Factoring<T>::ElementaryGuess::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(guess)
    ExportVarArrow(base)
    ExportVarArrow(weight)
    << "}";
}

template<typename T> auto
Factoring<T>::
Correction(const OffsetType & offset, int ord) const -> ScalarType {
    /*std::cout
    ExportArrayArrow(guesses)
    ExportVarArrow(offset)
    << "\n";*/
//    if(true || pointChoice!=FactoringPointChoice::Both){
        return std::accumulate(guesses.begin(),guesses.end(),0.,
                               [&](ScalarType a, const ElementaryGuess & b)->ScalarType{
                                   return a + b.Correction(offset,ord);});
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

template<typename T> bool
Factoring<T>::NeedsRecompute(IndexCRef index) const {
	if(method == FactoringMethod::None) return false;
	assert(!factoringRegion.empty());
	return factoringRegion(index);
}

template<typename T> bool
Factoring<T>::SetIndexStatic(IndexCRef index){
	assert(method==FactoringMethod::Static);
	assert(!factoringRegion.empty());
	const DiscreteType linearIndex = factoringRegion.Convert(index);
	if(!factoringRegion[linearIndex]) return false;
	
	assert(false);
	//TODO : complete here
	ExceptionMacro("Factoring error : Static factoring is unsupported for now")
}


template<typename T> bool // Returns wether dynamic factoring is active at this point
Factoring<T>::
SetIndexDynamic(IndexCRef index, const DiscreteFlowType & flow) {
	assert(method==FactoringMethod::Dynamic);
    assert(!factoringRegion.empty());
    const DiscreteType linearIndex = factoringRegion.Convert(index);
    if(!factoringRegion[linearIndex]) return false;

//    std::cout ExportVarArrow(index) << "\n";
    const FullIndexCRef full{index,linearIndex};
    MakeFactor(full, flow);
    return MakeGuess(full);
}

template<typename T> void
Factoring<T>::
MakeFactor(FullIndexCRef updated, const DiscreteFlowType & flow){
    if(factoringDone[updated.linear]) return; // Work already done
    factoringDone[updated.linear]=true;
    
    currentNeighbors.clear();
    currentNeighbors2.clear();
    
    
    const ScalarType wSum = std::accumulate(flow.begin(), flow.end(), 0.,
         [](ScalarType a, const DiscreteFlowElement & b)->ScalarType{return a+b.weight;});
    assert(wSum>0);
    
    for(const auto & [offset,weight] : flow){
        IndexType neigh = updated.index+IndexDiff::CastCoordinates(offset);
        assert(pFM!=nullptr);
        [[maybe_unused]] const auto transform = pFM->dom.Periodize(neigh,updated.index);
        assert(transform.IsValid());
        const DiscreteType linearNeigh = factoringRegion.Convert(neigh);
        
        if(factoringKeypoints.find(linearNeigh)!=factoringKeypoints.end()){
            currentNeighbors.push_back({offset,weight});
        } else {
            const auto [neigh_begin,neigh_end] = factoringNeighbors.equal_range(linearNeigh);
            for(auto it = neigh_begin; it!=neigh_end; ++it){
                const auto & [offset2,weight2] = it->second;
                currentNeighbors.push_back({offset+offset2, weight*weight2});
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
    for(const auto & [offset,weight] : currentNeighbors){
        if(currentNeighbors2.empty() || currentNeighbors2.back().first != offset){
			currentNeighbors2.push_back({offset,weight});
        } else {
            currentNeighbors2.back().second+=weight;
        }
    }
    
    const auto it = factoringNeighbors.upper_bound(updated.linear);
    for(const auto & data : currentNeighbors2){
        factoringNeighbors.insert(it,{updated.linear,data});}
}

template<typename T> bool
Factoring<T>::
MakeGuess(FullIndexCRef updated) {
    assert(factoringRegion[updated.linear]);
    assert(factoringDone[updated.linear]);
    guesses.clear();
    const auto rg = factoringNeighbors.equal_range(updated.linear);
    if(rg.first==rg.second) return false;
    
    if(pointChoice==FactoringPointChoice::Current || pointChoice==FactoringPointChoice::Both){
        const DistanceGuess guess = pFM->stencilData.GetGuess(updated.index);
        for(auto it = rg.first; it!=rg.second; ++it){
            const auto & [offset,weight] = it->second;
            const DomainTransformType transform{}; // Identity transform
            guesses.push_back({guess,transform,offset,
                weight *(pointChoice==FactoringPointChoice::Both ? 0.5 : 1.)});
        }
    }
    if(pointChoice==FactoringPointChoice::Key || pointChoice==FactoringPointChoice::Both){
        for(auto it = rg.first; it!=rg.second; ++it){
            const auto & [offset,weight] = it->second;
            IndexType neigh = updated.index+IndexDiff::CastCoordinates(offset);
            const DomainTransformType transform = pFM->dom.Periodize(neigh,updated.index);
            const auto keyIt = factoringKeypoints.find(factoringDone.Convert(neigh));
            assert(keyIt!=factoringKeypoints.end());
            guesses.push_back({keyIt->second,transform,offset,
                weight *(pointChoice==FactoringPointChoice::Both ? 0.5 : 1.)});
        }
    }
    return true;
}

template<typename T> bool
Factoring<T>::
Setup(HFMI * that){
    auto & io = that->io;
	
	// Get the factoring method
	
	method=enumFromString<FactoringMethod>
	(io.GetString("factoringMethod",enumToRealString(method)));
	if(0>(int)method){ExceptionMacro("Dynamic factoring error: unrecognized factoringMethod.");}
	if(method==FactoringMethod::None) {return false;}
	
	// Get the factoring point choice
	pointChoice = enumFromString<FactoringPointChoice>
	(io.GetString("factoringPointChoice",enumToRealString(pointChoice)));
	if(0>(int)pointChoice){
		ExceptionMacro("Dynamic factoring error: unrecognized factoringPointChoice.");}
	
	// Setup arrays

	factoringRegion.dims = pFM->stencilData.dims;
	factoringRegion.resize(factoringRegion.dims.Product(),false);
	if(method==FactoringMethod::Dynamic){
		factoringDone.dims = factoringRegion.dims;
		factoringDone.resize(factoringRegion.size(),false);
	}

	// Get the factoring points

	SetupCenters(that);
	SetupRegion(that);
	
	// Export some requested data
	if(io.template Get<ScalarType>("exportFactoringCenters",0.,2)){
		const auto & param = that->pFM->stencilData.Param();
		std::vector<PointType> pts;
		for(const auto & [p,d] : factoringCenters) {
			pts.push_back(param.ReDim(p));}
		io.SetVector("factoringCenters",pts);
	}

	if(io.template Get<ScalarType>("exportFactoringRegion",0.,2)!=0){
		io.SetArray("factoringRegion",factoringRegion.template Cast<ScalarType>());}

    return true;
}

template<typename T> void Factoring<T>::
SetupDijkstra(){
	if(!edges.empty()) return;
	// Setup the edges used in Dijkstra's method, coarsely approximating euclidean distance
	Array<int, Dimension> dummy; dummy.dims = IndexType::Constant(3);
	int nEdges = dummy.dims.Product();
	for(int i=0; i<nEdges; ++i){
		const IndexDiff offset = dummy.Convert(i)-IndexType::Constant(1);
		const ScalarType norm = VectorType::CastCoordinates(offset).Norm();
		if(norm>0){edges.push_back({offset,norm});}
	}
}


template<typename T> void Factoring<T>::
SetupRegion(HFMI*that){
	auto & io = that->io;
	pFM = that->pFM.get();
	const auto & dom = pFM->dom;

	// Get the factoring region map, if user provides
	
	if(method==FactoringMethod::Dynamic || !io.HasField("factoringRegion")){
		factoringRadius = io.template Get<ScalarType>("factoringRadius",factoringRadius);}
	
	if(io.HasField("factoringRegion")){
		const auto pRegion = that->template GetIntegralField<bool>("factoringRegion");
		pRegion->CheckDims(factoringRegion.dims);
		for(int i=0; i<factoringRegion.size(); ++i){
			factoringRegion[i] = (*pRegion)(factoringRegion.Convert(i));}
		return;
	}
	
	// Build the factoring region, based on the factoring radius and known centers
	const auto & arr = factoringRegion;
	
	// ------ Setup dynamic factoring region,
	// by running a Dijkstra, approximating the euclidean distance ------

	SetupDijkstra();
	
	std::priority_queue<std::pair<ScalarType,DiscreteType> > queue;
	for(const auto & [p,d] : factoringCenters){
		queue.push({0.,arr.Convert(dom.IndexFromPoint(p))});}
	
	while(!queue.empty()){
		const auto [value,current] = queue.top();
		queue.pop();
		if(factoringRegion[current]) continue;
		factoringRegion[current]=true;
		for(const auto & [offset,cost] : edges){
			const ScalarType neighVal = value+cost;
			if(neighVal>factoringRadius) continue;
			const IndexType index = arr.Convert(current);
			IndexType neigh = index+offset;
			if(dom.Periodize(neigh,index).IsValid()){
				queue.push({-neighVal,arr.Convert(neigh)});}
		}
	}
	
}


template<typename T> void Factoring<T>::
SetupCenters(HFMI * that) {
	
	IO & io = that->io;
	const auto & stencilData = pFM->stencilData;
	const auto & param = stencilData.Param();
	const auto & dom = pFM->dom;
	Array<bool, Dimension> & done = factoringRegion;
	assert(!done.empty());


	// --- import the used provided centers ---
	auto pushCenter = [this,&stencilData,&dom,&done](PointType q){
		if(!dom.PeriodizeNoBase(q).IsValid()){
			ExceptionMacro("Dynamic factoring error : input point " << q << " is out of range");}
		if(method==FactoringMethod::Dynamic){
			for(const auto [index,weight] : dom.Neighbors(q)){
				if(weight>1e-6){
					factoringKeypoints.insert({done.Convert(index),factoringCenters.size()});}
			}
		}
		factoringCenters.push_back({q,stencilData.GetGuess(q)});
	};
	
	for(const PointType & p : io.template GetVector<PointType>("seeds")) {
		pushCenter(param.ADim(p));}
	
	if(method!=FactoringMethod::Dynamic) return;
	if(io.HasField("factoringKeypoints")) {
		for(const PointType & p : io.template GetVector<PointType>("factoringKeypoints")){
			pushCenter(param.ADim(p));}
		return;
	}

    // In the case of dynamic factoring, without user provided keypoints,
	// we extract keypoints from wall data, trying to get the corners.
	// This relies on a Dijkstra method.
	if(!io.HasField("walls")) return;

	SetupDijkstra();

    const auto pWalls = that->template GetIntegralField<bool>("walls");
    pWalls->CheckDims(done.dims);
    assert(!edges.empty());
    
    const ScalarType factoringWallExtractionRadius = io.template Get<ScalarType>("factoringWallExtractionRadius",factoringRadius,2);
	
    // Get the wall boundary, which are the candidate keypoints
    std::set<DiscreteType> wallBoundary;
    for(DiscreteType i=0; i<done.size(); ++i){
        const IndexType index = done.Convert(i);
        if(!(*pWalls)(index)) continue;
        for(const auto & [offset,cost] : edges){
            IndexType neigh = index +offset;
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
        const auto [value,current,source] = queue.top();
        queue.pop();
        if(done[current]) continue;
        done[current]=true;
        if(source!=current){
            const auto it=score.find(source);
            if(it!=score.end()) {++(it->second);}
            else {score.insert({source,1});}
        }
        for(const auto & [offset,cost] : edges){
            const ScalarType neighVal = value+cost;
            if(neighVal>factoringWallExtractionRadius) continue;
            const IndexType index = done.Convert(current);
            IndexType neigh = index+offset;
            if(dom.Periodize(neigh,index).IsValid() && !(*pWalls)(neigh)){
                queue.push({-neighVal,done.Convert(neigh),source});}
        }
    }
	
    const ScalarType factoringWallExtractionThreshold = io.template Get<ScalarType>("factoringWallExtractionThreshold", factoringRadius*sqrt(factoringRadius), 2);
    for(const auto & [index,value] : score){
		if(value>factoringWallExtractionThreshold && // Singular point of the walls,
		   !OnBoundary(done.Convert(index), dom) ){// not on domain boundary
			pushCenter(dom.PointFromIndex(index));
		}
	}
	
    std::fill(done.begin(),done.end(),false); // reset this array for future use
}

template<typename T> bool
Factoring<T>::
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
