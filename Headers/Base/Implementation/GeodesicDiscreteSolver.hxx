// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef GeodesicDiscreteSolver_hxx
#define GeodesicDiscreteSolver_hxx

template<typename Traits> void
GeodesicDiscreteSolver<Traits>::Setup(HFMI * that) {
    assert(that!=nullptr);
    auto & io = that->io;
    geodesicStep = io.template Get<ScalarType>("geodesicStep",geodesicStep);
    weightThreshold = io.template Get<ScalarType>("geodesicWeightThreshold",weightThreshold);
    volumeBound = io.template Get<ScalarType>("geodesicVolumeBound",volumeBound);
}

template<typename Traits> std::vector<std::vector<typename Traits::PointType> >
GeodesicDiscreteSolver<Traits>::Run(HFMI * that, const std::vector<PointType> & tips) {
    assert(that!=nullptr);
    std::vector<std::vector<PointType> > result;
    for(const PointType & tip : tips){
        std::vector<PointType> geodesic;
        geodesic.push_back(tip);
        const bool failed = Run(geodesic);
        if(failed) Msg() << "Geodesic extraction failed for tip "
			<< that->stencil.Param().ReDim( tip ) << "\n";
        else if(nRestarts>nRestartsBeforeVolumeIncrease){
            Msg() << "Tip " << that->stencil.Param().ReDim( tip )
            << " yields " << nRestarts << " restarts, "
            << "geodesicVolumeBound increased to " << effectiveVolumeBound;}
        result.push_back(std::move(geodesic));
    }
	
    return result;
}

// ------------------------------------------
// -- Approach based on Stencil structure ---
// ------------------------------------------

template<typename Traits> bool GeodesicDiscreteSolver<Traits>::
Run(std::vector<PointType> & geodesic) const  {
    if(geodesic.empty()) ExceptionMacro("GeodesicDiscreteSolver error : no start point !");
        nRestarts = 0;
        effectiveVolumeBound = volumeBound;
        while(nRestarts<nMaxRestarts){
            const bool needsRestart = GeodesicDiscrete(geodesic);
            if(!needsRestart) break;
            ++nRestarts;
            if(nRestarts>=nRestartsBeforeVolumeIncrease){
                effectiveVolumeBound *= volumeIncreaseRatio;}
        }
    return nRestarts==nMaxRestarts;
}

template<typename Traits> bool GeodesicDiscreteSolver<Traits>::
GeodesicDiscrete(std::vector<PointType> & geodesic) const {
    const VectorType HalfVec = VectorType::Constant(0.5);
    const auto & fm = this->fm;
    const auto & dom = fm.dom;
    
    // (value,index) -> weight
    typedef std::map<std::pair<ScalarType,IndexType>,ScalarType> GeoType;
    GeoType geo;
    VectorType indexSum=VectorType::Constant(0);
    ScalarType weightSum=0;
    PointType indexAvg;
    typedef typename Traits::template BasisReduction<Dimension>::SymmetricMatrixType Sym;
    Sym covSum = Sym::Zero();
    
    auto UpdateSums = [&weightSum,&indexSum,&covSum](ScalarType weight, IndexType ind)->void {
        weightSum+=weight;
        const VectorType v = VectorType::FromOrigin(PointType::CastCoordinates(ind));
        indexSum += weight*v;
        covSum += weight*Sym::RankOneTensor(v);
    };
    	
    // (value, increment)
    std::vector<std::pair<ScalarType,VectorType> > steps;
    ScalarType stepNorm=0;
    PointType indexBase;
    
    { 	// Set a unit weight, split among the neighbors on the grid
        assert(!geodesic.empty());
        const PointType tip = geodesic.back();
        const IndexType indexRef = dom.IndexFromPoint(tip);
        const PointType pointRef = dom.PointFromIndex(indexRef); // Closest point on grid
        for(int i=0; i< (1<< Dimension); ++i){
            // TODO : find out if in wall, using visible offset
            ScalarType weight=1;
            IndexType index = indexRef;
            for(int j=0; j<Dimension; ++j){
                const ScalarType diff = tip[j]-pointRef[j];
                if( i & (1 << j) ){
                    index[j]+=diff>0 ? 1 : -1;
                    weight *= std::abs(diff);
                } else {
                    weight *= 1.-std::abs(diff);
                }
            } // for j
            assert(weight>=0);
            IndexType indexPer=index;
            if(!dom.Periodize(indexPer,indexRef).IsValid()) continue;
            geo[{fm.values(indexPer),index}] = weight;
            UpdateSums(weight,index);
        }
		if(weightSum==0)  {return false;}
        indexAvg = PointType::FromOrigin(indexSum/weightSum);
        indexBase= indexAvg;
    }
	
	// Weight sum should remain approximately one. Some loss can happen, but
	// vey small values are irrelevant. (Encountered with some cycling due to loss
	// of causality when using factoring without spreadSeeds)
	const ScalarType weightSumLowerBound = 1e-5;
	
    while(!geo.empty() && weightSum>weightSumLowerBound){
        // Extract and erase point with largest value.
        const auto it = --geo.end(); // point with largest value
		const auto [value,index] = it->first;
        const ScalarType weight = it->second;
        geo.erase(it);
				
        if(weight==0) continue;
        UpdateSums(-weight,index);
		
        // Insert children points.
        IndexType perIndex = index;
        const auto transform = dom.Periodize(perIndex,perIndex);
        assert(transform.IsValid());
        DiscreteFlowType flow;
        fm.Recompute(perIndex, flow);
		
        ScalarType wSum=0;
        for(const auto & offsetWeight : flow){
            wSum+=offsetWeight.weight;}
		if(wSum==0) continue;
        
        
        for(const auto & offsetWeight : flow){
            const ScalarType w= weight*offsetWeight.weight/wSum;
            IndexType ind = index;
            OffsetType offset = offsetWeight.offset;
            transform.PullVector(offset);
            ind+=IndexDiff::CastCoordinates(offset);
            IndexType perInd = ind;
            const auto transform = dom.Periodize(perInd,index);
            assert(transform.IsValid()); (void)transform;
            const ScalarType val = fm.values(perInd);
            const auto it = geo.find({val,ind});
            if(it==geo.end()){
                // Reject new point if too small weight, or too far.
				if(w<weightSum*weightThreshold)  continue;
                geo[{val,ind}] = w;
            } else {
                it->second+=w;}
            UpdateSums(w,ind);
        }
		        
        // Compute step
        if(weightSum<=0) {
            geodesic.push_back(indexAvg+HalfVec);
            return false;}
        const PointType oldIndexAvg = indexAvg;
        indexAvg = PointType::FromOrigin(indexSum/weightSum);
        const VectorType step = indexAvg-oldIndexAvg;
        steps.push_back({value,step});
        stepNorm+=step.Norm();
								 
        if(geo.size()==1){ // If geodesic concentrates at a seed or around a wall corner.
			if(flow.empty()) break; // At seed
            geodesic.push_back(indexAvg+HalfVec);
            indexBase=indexAvg;
            stepNorm=0;
            steps.clear();			
            continue;
        }
        
        // Perform step, if necessary
        if(stepNorm<geodesicStep && geo.size()>1) continue;
        const ScalarType valInc =  steps.front().first-value;
        if(valInc<=0) continue;
        const ScalarType delta = valInc/step.size();
        PointType p = indexBase;
        for(const auto & valStep : steps){
            p+=valStep.second * (valStep.first-value+delta)/(valInc+delta);}
        geodesic.push_back(p+HalfVec);
		
        indexBase = indexAvg;
        stepNorm=0;
        steps.clear();
        
        
        // If points are too much spread, restart
        const Sym variance = covSum/weightSum-Sym::RankOneTensor(indexSum/weightSum)+Sym::Identity();
        const ScalarType detVar = variance.Determinant();
        
        if(detVar>square(effectiveVolumeBound)){
            geodesic.push_back(indexAvg+HalfVec);
            return true;
        } /*else {
		   // Do some cleanup of faraway points with little weight ?
            // Seems to have little impact on results
		   
		   // Problematic
		   // Too strict criteria on large instances.
		   // Eliminates too much mass.
		   // Also seems to have unreasonable computation time.
            
            // Alternative pruning method tried, based on anisotropic spread.
            // Quite costly, and not very efficient.
		   
            if(Dimension>3) continue;
            const ScalarType meanWeight = weightSum/geo.size();
            const Sym invVar = variance.Inverse();
            for(auto & indVal : geo){
                const IndexType ind = indVal.first.second;
                ScalarType & w = indVal.second;
                if(w>=meanWeight/2) continue;
                const ScalarType sqNorm = invVar.SquaredNorm(PointType::CastCoordinates(ind)-indexAvg);
                if(detVar*sqNorm<=square(2.*effectiveVolumeBound)) continue;
                UpdateSums(-w,ind);
                w=0;
            }
            indexAvg=PointType::FromOrigin(indexSum/weightSum);
			
        }*/
    }
	
    return false;
}

#endif /* GeodesicDiscreteSolver_hxx */
