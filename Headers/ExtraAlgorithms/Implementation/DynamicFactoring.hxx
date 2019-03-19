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

template<typename T> void Factoring<T>::
ElementaryGuess::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(centerIndex)
    ExportVarArrow(base)
    << "}";
}

template<typename T> auto Factoring<T>::
Correction(const OffsetType & offset, int ord) const -> ScalarType {
	
	/*
	// ----- Debug with exact value for Poincare Half disk test ---
	// Scale is 0.1,
	const PointType seed(0.,1.5);
	auto Value = [&seed](const PointType & q){
		const VectorType v = q-seed;
		const VectorType w (q[0]-seed[0],q[1]+seed[1]);
		return log( (w.Norm()+v.Norm()) / (w.Norm()-v.Norm()) );
	};
	auto Gradient = [&seed](const PointType & q){
		const VectorType v = q-seed;
		const VectorType w (q[0]-seed[0],q[1]+seed[1]);
		const ScalarType wN = w.Norm(), vN = v.Norm();
		const VectorType wU = w/wN, vU = v/vN;
		return (wU+vU)/(wN+vN) - (wU-vU)/(wN-vN);
	};
	
	const ScalarType h = 0.1;
	const PointType  q = seed - h * data_static.base;
	const VectorType o = h*VectorType::CastCoordinates(offset);
	
	
	std::cout << "In correction "
	ExportVarArrow(q)
	ExportVarArrow(q+o)
	ExportVarArrow(Value(q))
	ExportVarArrow(Value(q+o))
	ExportVarArrow(Gradient(q))
	<< std::endl;
	
	
	switch(ord){
			//case 1: return deriv+(valueZero-Value(offset));
		case 1: return Gradient(q).ScalarProduct(o) + (Value(q)-Value(q+o));
		default:assert(false);
	}
	*/
	
	// ---------
	
	
	auto Correct = [&](const ElementaryGuess & data){
		switch(this->pointChoice){
			case FactoringPointChoice::Current:
				return CorrectionCurrent(offset, ord, data.base);
			case FactoringPointChoice::Key:
				return CorrectionKey(offset,ord, data);
			case FactoringPointChoice::Both:
				return 0.5*(CorrectionCurrent(offset, ord, data.base)
							+CorrectionKey(offset,ord, data));
		} // Switch point
	};
	
	if(method==FactoringMethod::Static && data_static.centerIndex != -1){
		return Correct(data_static);

		//DEBUG
		useExact=false;
		const ScalarType correction_inexact = Correct(data_static);
		useExact=true;
		const ScalarType correction_exact = Correct(data_static);
		const ScalarType correction_rdiff = (correction_inexact-correction_exact)/correction_exact;
		std::cout << "In correction "
		ExportVarArrow(correction_inexact)
		ExportVarArrow(correction_exact)
		ExportVarArrow(correction_rdiff)
		ExportVarArrow(data_static.base)
		ExportVarArrow(offset)
		ExportVarArrow(ord)
		<< std::endl << std::endl;
		return correction_exact;
		
		return Correct(data_static);
	} else if(method==FactoringMethod::Dynamic && !data_dynamic.guesses.empty() ){
		return std::accumulate(data_dynamic.guesses.begin(),data_dynamic.guesses.end(),0.,
				[&](ScalarType a, const std::pair<ScalarType,ElementaryGuess> & b)->ScalarType{
					return a + b.first*Correct(b.second);});
	} else {
		return 0.;
	}
}

// ---------------- DEBUG -----------------
template<typename T>
struct Factoring<T>::PoincareExact {
	const PointType seed = PointType(0.,1.5);
	const ScalarType h = 0.1;
	ScalarType Value(const PointType & q) const {
		const VectorType v = q-seed;
		const VectorType w (q[0]-seed[0],q[1]+seed[1]);
		return log( (w.Norm()+v.Norm()) / (w.Norm()-v.Norm()) );
	};
	VectorType Gradient(const PointType & q) const {
		const VectorType v = q-seed;
		const VectorType w (q[0]-seed[0],q[1]+seed[1]);
		const ScalarType wN = w.Norm(), vN = v.Norm();
		const VectorType wU = w/wN, vU = v/vN;
		return h*((wU+vU)/(wN+vN) - (wU-vU)/(wN-vN));
	};
	PointType Pos(const VectorType &base) const {
		return seed-h*base;
	}
};


template<typename T> auto Factoring<T>::
CorrectionCurrent(const OffsetType & offset,
				  int order,
				  const VectorType & base
) const -> ScalarType {
/*	if(order>=2){
	std::cout << "In Correction current"
	ExportVarArrow(order)
	ExportVarArrow(base)
	ExportVarArrow(offset)
	<< std::endl;
	}*/
	
	assert(pointChoice==FactoringPointChoice::Current || pointChoice==FactoringPointChoice::Both);
	const DistanceGuess & dist = this->currentGuess;
	auto Value = [&base,&dist](const OffsetType & off) -> ScalarType {
		const VectorType v = base-VectorType::CastCoordinates(off);
		return v.IsNull() ? 0. : dist.Norm(v);
	};
	const ScalarType valueZero = dist.Norm(base);
	const ScalarType deriv = -
	dist.Gradient(base).ScalarProduct(VectorType::CastCoordinates(offset));
	
	
	if(useExact){// DEBUG
		PoincareExact exact;
		auto Value = [&](VectorType off){
			return exact.Value(exact.Pos(base-off ));
		};
		const VectorType pulledOffset = VectorType::CastCoordinates(offset);
		const ScalarType deriv = exact.Gradient(exact.Pos(base)).ScalarProduct(pulledOffset);
		const ScalarType valueZero = Value(VectorType::Constant(0));
		switch(order){
			case 1: return deriv+ (valueZero - Value(pulledOffset));
			case 2: return deriv+(1.5*valueZero - 2.*Value(pulledOffset) + 0.5*Value(2*pulledOffset));
			case 3: return deriv
				+ ((11./6.)*valueZero-3.*Value(pulledOffset)+1.5*Value(2*pulledOffset)-(1./3.)*Value(3*pulledOffset));
				
		}
	} 	// End DEBUG
	
	
	
	switch(order){
		case 1: return deriv+(valueZero-Value(offset));
		case 2: return deriv+(1.5*valueZero - 2.*Value(offset) + 0.5*Value(2*offset));
		case 3: {
			return deriv
			+ ((11./6.)*valueZero-3.*Value(offset)+1.5*Value(2*offset)-(1./3.)*Value(3*offset));
			
			// Below is old, erroneous implem
			// According to further calculations, the spatial derivative of the metric should not be added, since we do not change the base point in the Value lambda function.
			
			
			IndexType neigh = currentIndex+IndexDiff::CastCoordinates(offset);
			const DomainTransformType transform = pFM->dom.Periodize(neigh, currentIndex);
			assert(transform.IsValid());
			VectorType pulledBase = base; transform.PullVector(pulledBase);
			const DistanceGuess neighGuess = pFM->stencilData.GetGuess(neigh);
			
			// Testing with a second neighbor
			IndexType neigh2 = currentIndex+2*IndexDiff::CastCoordinates(offset);
			const DistanceGuess neighGuess2 = pFM->stencilData.GetGuess(neigh2);
			
			const ScalarType diff = neighGuess.Norm(pulledBase) - valueZero;
			const ScalarType diff2 = (-0.5*neighGuess2.Norm(pulledBase)+2*neighGuess.Norm(pulledBase) - 1.5*valueZero);
			std::cout << "In Correction current - 3 : "
			ExportVarArrow(base)
			ExportVarArrow(offset)
			ExportVarArrow(deriv)
			ExportVarArrow(diff)
			ExportVarArrow(diff2)
			<< std::endl
			ExportVarArrow(dist)
			ExportVarArrow(neighGuess)
			ExportVarArrow(neighGuess2)
			<< std::endl;
			
			
			return deriv
			+ ((11./6.)*valueZero-3.*Value(offset)+1.5*Value(2*offset)-(1./3.)*Value(3*offset))
			+ (-0.5*neighGuess2.Norm(pulledBase)+2*neighGuess.Norm(pulledBase) - 1.5*valueZero);
			
			return deriv
			+ ((11./6.)*valueZero-3.*Value(offset)+1.5*Value(2*offset)-(1./3.)*Value(3*offset))
			+ (neighGuess.Norm(pulledBase) - valueZero);
		}
		default: assert(false); return -Traits::Infinity();
	}
}

template<typename T> auto Factoring<T>::
CorrectionKey(const OffsetType & offset, int order, const ElementaryGuess & guess)
const -> ScalarType {
	assert(0<=guess.centerIndex && guess.centerIndex<factoringCenters.size());
	const auto & dist = factoringCenters[guess.centerIndex].second;
	const auto & transform = guess.transform;

	VectorType pulledBase = guess.base; transform.PullVector(pulledBase);
	VectorType pulledOffset = VectorType::CastCoordinates(offset); transform.PullVector(pulledOffset);
	
	auto Value = [&pulledBase,&dist](const VectorType & off) -> ScalarType {
		const VectorType v=pulledBase-off;
		return v.IsNull() ? 0. : dist.Norm(v);
	};
	const ScalarType valueZero = dist.Norm(pulledBase);
	const ScalarType deriv = -dist.Gradient(pulledBase).ScalarProduct(pulledOffset);
	
	if(useExact){// DEBUG
		PoincareExact exact;
		auto Value = [&](VectorType off){
			return exact.Value(exact.Pos(guess.base-off ));
		};
		const ScalarType deriv = exact.Gradient(exact.Pos(guess.base)).ScalarProduct(pulledOffset);
		const ScalarType valueZero = Value(VectorType::Constant(0));
		switch(order){
			case 1: return deriv+ (valueZero - Value(pulledOffset));
			case 2: return deriv+(1.5*valueZero - 2.*Value(pulledOffset) + 0.5*Value(2*pulledOffset));
			case 3: return deriv
				+ ((11./6.)*valueZero-3.*Value(pulledOffset)+1.5*Value(2*pulledOffset)-(1./3.)*Value(3*pulledOffset));

		}
	} 	// End DEBUG

	
	switch(order){
		case 1:	return deriv+(valueZero-Value(pulledOffset));
		case 2: return deriv+(1.5*valueZero - 2.*Value(pulledOffset) + 0.5*Value(2*pulledOffset));
		case 3:
			/*
			std::cout << "In correction key "
			ExportVarArrow(pulledBase)
			ExportVarArrow(offset)
			ExportVarArrow(deriv)
			<< std::endl;
			*/
			return deriv
			+ ((11./6.)*valueZero-3.*Value(pulledOffset)+1.5*Value(2*pulledOffset)-(1./3.)*Value(3*pulledOffset));
		default: assert(false); return -Traits::Infinity();
	}
}


template<typename T> bool Factoring<T>::
NeedsRecompute(IndexCRef index) {
	switch(method){
		case FactoringMethod::None: return false;
			
		case FactoringMethod::Dynamic: data_dynamic.guesses.clear();
		case FactoringMethod::Static:
		default:
			assert(!factoringRegion.empty());
			return factoringRegion(index);
	}
}

template<typename T> bool Factoring<T>::
SetIndexStatic(IndexCRef index){
	assert(method==FactoringMethod::Static);
	assert(!factoringRegion.empty());
	const DiscreteType linearIndex = factoringRegion.Convert(index);
	if(factoringRegion[linearIndex]) {
		
		if(pointChoice!=FactoringPointChoice::Key){
			currentIndex = index;
			currentGuess = pFM->stencilData.GetGuess(index);
		}
		
		// Find center within factoring radius
		assert(!factoringCenters.empty());
		const PointType p = pFM->dom.PointFromIndex(index);
		for(int i=0; i<factoringCenters.size(); ++i){
			const PointType & q = factoringCenters[i].first;
			const VectorType v = q-p;
			if(v.SquaredNorm() < square(factoringRadius)){
				data_static.centerIndex = i;
				data_static.base = v;
				return true;
			}
		}
	}
	
	data_static.centerIndex=-1;
	return false;
}


template<typename T> bool // Returns wether dynamic factoring is active at this point
Factoring<T>::
SetIndexDynamic(IndexCRef index, const DiscreteFlowType & flow) {
	assert(method==FactoringMethod::Dynamic);
    assert(!factoringRegion.empty());
    const DiscreteType linearIndex = factoringRegion.Convert(index);
    if(!factoringRegion[linearIndex]) return false;

	if(pointChoice!=FactoringPointChoice::Key){
		currentIndex = index;
		currentGuess = pFM->stencilData.GetGuess(index);
	}
	
    const FullIndexCRef full{index,linearIndex};
    MakeFactor(full, flow);
    return MakeGuess(full);
}

template<typename T> void
Factoring<T>::
MakeFactor(FullIndexCRef updated, const DiscreteFlowType & flow){
	assert(method==FactoringMethod::Dynamic);
	auto & data = data_dynamic;
    if(data.factoringDone[updated.linear]) return; // Work already done
    data.factoringDone[updated.linear]=true;
    
    data.currentNeighbors.clear();
    data.currentNeighbors2.clear();
	
    const ScalarType wSum = std::accumulate(flow.begin(), flow.end(), 0.,
         [](ScalarType a, const DiscreteFlowElement & b)->ScalarType{return a+b.weight;});
    assert(wSum>0);
    
    for(const auto & [offset,w] : flow){
		const ScalarType weight = w/wSum;
        IndexType neigh = updated.index+IndexDiff::CastCoordinates(offset);
        assert(pFM!=nullptr);
        const auto transform = pFM->dom.Periodize(neigh,updated.index);
        assert(transform.IsValid());
        const DiscreteType linearNeigh = factoringRegion.Convert(neigh);
        
        if(data.factoringKeypoints.find(linearNeigh)!=data.factoringKeypoints.end()){
            data.currentNeighbors.push_back({offset,weight});
        } else {
            const auto [neigh_begin,neigh_end] = data.factoringNeighbors.equal_range(linearNeigh);
            for(auto it = neigh_begin; it!=neigh_end; ++it){
                auto [offset2,weight2] = it->second;
				transform.PullVector(offset2);
                data.currentNeighbors.push_back(
				{offset+offset2, weight*weight2});
                // TODO : remove some of these, based on offset radius ?
            }
        }
    }
    
	assert(1.001 > std::accumulate(data.currentNeighbors.begin(),data.currentNeighbors.end(),0.,
								   [](ScalarType a, const std::pair<OffsetType,ScalarType> & b) -> ScalarType {return a+b.second;}) );
	
    // Merge duplicates and insert
    std::sort(data.currentNeighbors.begin(),data.currentNeighbors.end());
    for(const auto & [offset,weight] : data.currentNeighbors){
        if(data.currentNeighbors2.empty() || data.currentNeighbors2.back().first != offset){
			data.currentNeighbors2.push_back({offset,weight});
        } else {
            data.currentNeighbors2.back().second+=weight;
        }
    }
    
    const auto it = data.factoringNeighbors.upper_bound(updated.linear);
    for(const auto & neigh : data.currentNeighbors2){
        data.factoringNeighbors.insert(it,{updated.linear,neigh});}
}

template<typename T> bool
Factoring<T>::
MakeGuess(FullIndexCRef updated) {
	
	assert(method==FactoringMethod::Dynamic);

	const auto & data = data_dynamic;
	auto & guesses = data_dynamic.guesses;
    assert(factoringRegion[updated.linear]);
    assert(data.factoringDone[updated.linear]);
    const auto rg = data.factoringNeighbors.equal_range(updated.linear);
    if(rg.first==rg.second) return false;
	
	guesses.clear();
	for(auto it = rg.first; it!=rg.second; ++it){
		const auto & [offset,weight] = it->second;
		IndexType neigh = updated.index+IndexDiff::CastCoordinates(offset);
		const DomainTransformType transform = pFM->dom.Periodize(neigh,updated.index);
		const auto factIt = data.factoringKeypoints.find(data.factoringDone.Convert(neigh));
		assert(factIt!=data.factoringKeypoints.end());
		const int centerIndex = factIt->second;
		VectorType subOffset = factoringCenters[centerIndex].first - pFM->dom.PointFromIndex(neigh);
		transform.PullVector(subOffset);
		guesses.push_back({weight, ElementaryGuess{
			centerIndex,
			VectorType::CastCoordinates(offset)+subOffset,
			transform
			} });
	}
	
	/*
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
            const auto keyIt = data.factoringKeypoints.find(data.factoringDone.Convert(neigh));
            assert(keyIt!=data.factoringKeypoints.end());
            guesses.push_back({factoringCenters[keyIt->second].second,transform,offset,
                weight *(pointChoice==FactoringPointChoice::Both ? 0.5 : 1.)});
        }
    }
	 */
    return true;
}

template<typename T> bool
Factoring<T>::
Setup(HFMI * that){
    auto & io = that->io;
	pFM = that->pFM.get();
	// Get the factoring method
	
	method=enumFromString<FactoringMethod>
	(io.GetString("factoringMethod",enumToRealString(method)));
	if(0>(int)method){ExceptionMacro("Dynamic factoring error: unrecognized factoringMethod.");}
	if(method==FactoringMethod::None) {return false;}
	
	// Get the factoring point choice
	if(pFM->order==3) pointChoice=FactoringPointChoice::Both; // Appropriate default for third order
	pointChoice = enumFromString<FactoringPointChoice>
	(io.GetString("factoringPointChoice",enumToRealString(pointChoice)));
	if(0>(int)pointChoice){
		ExceptionMacro("Dynamic factoring error: unrecognized factoringPointChoice.");}
	
	// Setup arrays

	factoringRegion.dims = pFM->stencilData.dims;
	factoringRegion.resize(factoringRegion.dims.Product(),false);
	if(method==FactoringMethod::Dynamic){
		data_dynamic.factoringDone.dims = factoringRegion.dims;
		data_dynamic.factoringDone.resize(factoringRegion.size(),false);
	}

	// Get the factoring points

	SetupCenters(that);
	SetupRegion(that);

	// Export some requested data
	if(io.template Get<ScalarType>("exportFactoringCenters",0.,2)!=0){
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
		const auto [minusValue,current] = queue.top();
		queue.pop();
		if(factoringRegion[current]) continue;
		factoringRegion[current]=true;
		for(const auto & [offset,cost] : edges){
			const ScalarType neighVal = -minusValue+cost;
			if(neighVal>factoringRadius) continue;
			const IndexType index = arr.Convert(current);
			IndexType neigh = index+IndexDiff::CastCoordinates(offset);
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
					data_dynamic.factoringKeypoints.insert(
					{done.Convert(index),factoringCenters.size()});}
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
    for(const auto & [linearIndex,value] : score){
		const IndexType & index = done.Convert(linearIndex);
		if(value>factoringWallExtractionThreshold && // Singular point of the walls,
		   !OnBoundary(index, dom) ){// not on domain boundary
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
