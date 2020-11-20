// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef FirstVariation_h
#define FirstVariation_h

// Forward and backward differentiation of the value with respect to variations of the cost function
// (which is the inverse of the speed function), and of the seeds values.
// TODO : make this exact as well in the case of second order differences/time varying speed field.

template<typename T> struct FirstVariation :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    Redeclare12Types(HFM,IndexType,ScalarType,IndexCRef,HFMI,DifferenceType,
					 ActiveNeighFlagType,DiscreteType,Traits,RecomputeType,
					 DiscreteFlowType,PointType,IndexDiff)
    Redeclare1Constant(HFM,Dimension)

    template<typename E, size_t n> using Array = typename HFM::template Array<E,n>;
    typedef typename std::conditional<HFM::policy==SSP::Share, typename HFM::MultiplierType, ScalarType>::type MultType;
    
    ScalarType * pCurrentTime=nullptr;
    struct NeighborType {IndexType index; ScalarType weight;};
    typedef CappedVector<NeighborType, HFM::StencilType::nActiveNeigh> ActiveNeighborsType;
    MultType ValueVariation(IndexCRef, ActiveNeighborsType &) const; // Elementary differentiation
    std::vector<std::pair<IndexType,ScalarType> >
    BackwardVariation(const std::vector<std::pair<IndexType,ScalarType> > &, Array<MultType, Dimension> &) const;
    void ForwardVariation(const Array<MultType, Dimension+1> &, Array<ScalarType, Dimension+1> &) const;
	
	virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*_pFM) override {pFM=_pFM; return true;}
protected:
    const HFM * pFM;
	static const size_t multSize0 = DifferenceType::multSize >=0 ? DifferenceType::multSize : 0;
    template<size_t VMultSize=multSize0, typename Dummy=void> struct _DiffHelper;
    typedef _DiffHelper<> DiffHelper;
    typedef typename Array<ScalarType, Dimension+1>::IndexType DeepIndexType;
    DeepIndexType DeepIndex(IndexType index, DiscreteType k) const {
        DeepIndexType result;
        for(int i=0; i<Dimension; ++i) {result[i]=index[i];}
        result[Dimension]=k;
        return result;
    };
    template<bool b=( HFM::policy==SSP::Lag2 || HFM::policy==SSP::Lag3 ), typename Dummy=void> struct ValueVariationHelper;
};

// --------- Setup ---------
template<typename T> void FirstVariation<T>::Setup(HFMI*that){
	auto & io = that->io;
	if(io.HasField("inspectSensitivity") || io.HasField("costVariation") || io.HasField("seedValueVariation")){
		if(that->pFM->factoring.seedRadius!=0)
			WarnMsg() << "First variation warning : spread seeds are not supported.\n";
	}
}

// --------- Export -------

template<typename T> void FirstVariation<T>::Finally(HFMI*that){
    auto & io = that->io;
    if(io.HasField("inspectSensitivity")){
        const std::vector<PointType> inds = io.template GetVector<PointType>("inspectSensitivity");
        
        std::vector<ScalarType> lengths;
        if(io.HasField("inspectSensitivityLengths")){
            lengths=io.template GetVector<ScalarType>("inspectSensitivityLengths");
            if(std::accumulate(lengths.begin(),lengths.end(),0.)!=inds.size())
                ExceptionMacro("Error : inconsistent inspectSensitivityLengths");
        } else lengths.resize(inds.size(),1);
        
        std::vector<ScalarType> weights;
        if(io.HasField("inspectSensitivityWeights")) {
            weights=io.template GetVector<ScalarType>("inspectSensitivityWeights");
            if(weights.size()!=inds.size()) ExceptionMacro("Error: Inconsistent inspectSensitivityWeights size");
        } else weights.resize(inds.size(),1);
        
        std::vector<std::pair<IndexType, ScalarType> > iw;
        auto indIt=inds.begin(); auto wIt=weights.begin();
        Array<MultType,Dimension> sensitivity;
        for(int i=0; i<lengths.size(); ++i){
            for(int k=0; k<lengths[i]; ++k, ++indIt, ++wIt){
                PointType p = that->stencil.Param().ADim(*indIt);
                if(!that->pFM->dom.Periodize(p,p).IsValid()) {
                    ExceptionMacro("Error : inspectSensitivity data points are out of range");}
                iw.push_back({that->pFM->dom.IndexFromPoint(p),*wIt});
            }
            sensitivity.clear();
            auto sensitiveSeedIndices = BackwardVariation(iw,sensitivity);
            io.SetArray("costSensitivity_"+std::to_string(i),sensitivity);
//            MultArrayIO<>::Set(that,"costSensitivity_"+std::to_string(i),sensitivity);
            
            
            std::vector<std::pair<PointType,ScalarType> > sensitiveSeeds;
            sensitiveSeeds.reserve(sensitiveSeedIndices.size());
            for(const auto & iS : sensitiveSeedIndices){
                sensitiveSeeds.push_back({that->stencil.Param().ReDim(that->pFM->dom.PointFromIndex(iS.first)),iS.second});}
            io.SetVector("seedSensitivity_"+std::to_string(i),sensitiveSeeds);
            
            iw.clear();
        }
    }
    if(io.HasField("costVariation") || io.HasField("seedValueVariation")){
        Array<MultType,Dimension+1> fwdVar;
        if(io.HasField("costVariation")) fwdVar = io.template GetArray<MultType,Dimension+1>("costVariation");
        Array<ScalarType,Dimension+1> valueVariations;
        if(io.HasField("seedValueVariation")){
            const auto seeds = io.template GetVector<PointType>("seeds");
            const auto seedVariations = io.template GetArray<ScalarType,2>("seedValueVariation",IO::ArrayOrdering::RowMajor);
            
            const DiscreteType nVar = seedVariations.dims[0];
            valueVariations.dims=DeepIndex(that->stencil.dims, nVar);
            valueVariations.resize(valueVariations.dims.Product(),0);
            
            if(seeds.size()!=seedVariations.dims[1])
                ExceptionMacro("Error : Inconsistent size of seedValueVariation");
            for(int i=0; i<seeds.size(); ++i){
                const IndexType index = that->pFM->dom.IndexFromPoint(that->stencil.Param().ADim(seeds[i]));
                for(int k=0; k<nVar; ++k){valueVariations(DeepIndex(index, k)) = seedVariations({k,i});}
            }
        }
        ForwardVariation(fwdVar,valueVariations);
        io.SetArray("valueVariation", valueVariations);
    }
};

// ------------------- Distance gradient --------------------

template<typename TTraits> template<size_t VMultSize, typename Dummy>
struct FirstVariation<TTraits>::_DiffHelper {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare3Types(HFM,ScalarType,DifferenceType,MultiplierType);
    static MultiplierType NullMult(){return MultiplierType::Constant(0);}
    static ScalarType & Elem(MultiplierType & m, const DifferenceType & diff){return m[diff.multIndex];}
    static MultiplierType Times(ScalarType a, MultiplierType m){
        MultiplierType r; for(int i=0; i<VMultSize; ++i) r[i]=a*m[i]; return r;}
    static MultiplierType Times(const MultiplierType & m0, const MultiplierType & m1){
        MultiplierType r; for(int i=0; i<VMultSize; ++i) r[i]=m0[i]*m1[i]; return r;}
    static ScalarType Scal(const MultiplierType & m0, const MultiplierType & m1){
        ScalarType r=0; for(int i=0; i<VMultSize; ++i) r+=m0[i]*m1[i]; return r;}
};


template<typename TTraits> template<typename Dummy>
struct FirstVariation<TTraits>::_DiffHelper<1,Dummy> {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare2Types(HFM,ScalarType,DifferenceType);
    static ScalarType NullMult(){return 0;}
    static ScalarType & Elem(ScalarType & m, const DifferenceType & diff){return m;}
    static ScalarType Times(ScalarType a, ScalarType m){return a*m;}
    static ScalarType Scal(ScalarType m0, ScalarType m1){return m0*m1;}
};

template<typename TTraits> template<typename Dummy>
struct FirstVariation<TTraits>::_DiffHelper<0,Dummy> {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare3Types(HFM,ScalarType,DifferenceType,MultiplierType);
    static ScalarType NullMult(){return 0;}
    static ScalarType & Elem(ScalarType & m, const DifferenceType & diff){return m;}
    static ScalarType Times(ScalarType a, ScalarType m){return a*m;}
    static ScalarType Times(MultiplierType a, ScalarType m){return m;}
    static ScalarType Scal(ScalarType m0, ScalarType m1){return m0*m1;}
};

// ----------- Differentiation at a single point -----------



template<typename T> template<typename Dummy> struct
FirstVariation<T>::ValueVariationHelper<true,Dummy> {
    const HFM * pFM;
    MultType operator()(IndexCRef index, ActiveNeighborsType & neighbors) const {
        assert(neighbors.empty());
        DiscreteFlowType flow;
        const RecomputeType rec = pFM->Recompute(index, flow);
		
		typedef typename DiscreteFlowType::value_type FlowElementType;
		const ScalarType weightSum = std::accumulate(flow.begin(), flow.end(), 0.,[](ScalarType a, const FlowElementType & b){return a+b.weight;});
		
        ScalarType valueDiff = rec.value; // equivalently pFM->values(index)
        // valueDiff accounts for f(xMin) in the semi-Lagrangian optimization inf_x f(x) + u(x)
        // over the neighborhood boundary
        
        for(const auto & flowElem : flow){
            const ScalarType weight = flowElem.weight / weightSum;
            IndexType neigh = index+IndexDiff::CastCoordinates(flowElem.offset);
            [[maybe_unused]] const auto transform = pFM->dom.Periodize(neigh, index);
            assert(transform.IsValid());
            neighbors.push_back({neigh,weight});
            valueDiff -= weight*pFM->values(neigh);
        }
        
        assert(valueDiff>=0);
        return MultType{valueDiff};
    }
};

template<typename T> template<typename Dummy> struct
FirstVariation<T>::ValueVariationHelper<false,Dummy> {
    const HFM * pFM;
    MultType operator()(IndexCRef index, ActiveNeighborsType & neighbors) const {
        
        assert(neighbors.empty());
        const auto & values = pFM->values;
        const ScalarType value = values(index);
        
        MultType var = DiffHelper::NullMult();
        ScalarType weightSum=0;
        const auto & data = pFM->stencilData.RecomputeData(index);
        const ActiveNeighFlagType active = pFM->activeNeighs(index);
        
        auto func = [&,this](DiscreteType s, const DifferenceType & diff){
            IndexType neighIndex=index;
            neighIndex+=s*IndexDiff::CastCoordinates(diff.offset);
            const auto transform = pFM->dom.Periodize(neighIndex,index);
            assert(transform.IsValid()); (void)transform;
            const ScalarType neighValue = values(neighIndex);
            const ScalarType valueDiff = std::max(0., value-neighValue);
            const ScalarType w = diff.Weight(data.mult)*valueDiff;
            weightSum+=w;
            neighbors.push_back({neighIndex,w});
            DiffHelper::Elem(var,diff)+=w*valueDiff;
        };
        
        const int iMax = HFM::StencilType::GetIMax(active);
        int iNeigh = iMax*HFM::StencilType::nSingleNeigh;
        const auto & forward = data.stencil.forward[iMax];
        const auto & symmetric = data.stencil.symmetric[iMax];
        
        for(const auto & diff : forward){
            if(active[iNeigh]){
                func( 1,diff);}
            ++iNeigh;
        }
        for(const auto & diff : symmetric){
            if(active[iNeigh]){
                assert(!active[iNeigh+1]);
                func( 1,diff);}
            if(active[iNeigh+1]){
                assert(!active[iNeigh]);
                func(-1,diff);}
            iNeigh+=2;
        }
        
        
        if(weightSum==0) return var;
        const ScalarType invWeight = 1./weightSum;
        for(auto & neigh : neighbors){neigh.weight*=invWeight;}
        return DiffHelper::Times(invWeight,DiffHelper::Times(data.mult,var));
    }
};


template<typename T> auto FirstVariation<T>::
ValueVariation(IndexCRef index, ActiveNeighborsType & neighbors) const ->MultType {
    return ValueVariationHelper<>{pFM}(index,neighbors);
}

/*
template<typename T> auto FirstVariation<T>::
ValueVariation(IndexCRef index, ActiveNeighborsType & neighbors) const ->MultType {
    assert(neighbors.empty());
    DiscreteFlowType flow;
    const RecomputeType rec = pFM->Recompute(index, flow);
    
    typedef typename DiscreteFlowType::value_type FlowElementType;
    assert(std::fabs(std::accumulate(flow.begin(), flow.end(), 0.,[](ScalarType a, const FlowElementType & b){return a+b.weight;})-1.)<0.001);

    ScalarType valueDiff = rec.value; // equivalently pFM->values(index)
    // valueDiff accounts for f(xMin) in the semi-Lagrangian optimization inf f(xMin) + u(x)
    // over the neighborhood boundary
    
    for(const auto & flowElem : flow){
        const ScalarType weight = flowElem.weight;
        IndexType neigh = index+IndexDiff::CastCoordinates(flowElem.offset);
        const auto transform = pFM->dom.Periodize(neigh, index);
        assert(transform.IsValid());
        neighbors.push_back({neigh,weight});
        valueDiff -= weight*pFM->values(neigh);
    }
    
    assert(valueDiff>=0);
    return MultType{valueDiff};
}*/

/*
template<typename T> auto FirstVariation<T>::
ValueVariation(IndexCRef index, ActiveNeighborsType & neighbors) const ->MultType {
    assert(neighbors.empty());
    const auto & values = pFM->values;
    const ScalarType value = values(index);
    
    MultType var = DiffHelper::NullMult();
    ScalarType weightSum=0;
    const auto & data = pFM->stencilData.RecomputeData(index);
    const ActiveNeighFlagType active = pFM->activeNeighs(index);
    
    auto func = [&,this](DiscreteType s, const DifferenceType & diff){
        IndexType neighIndex=index;
        neighIndex+=s*IndexDiff::CastCoordinates(diff.offset);
        const auto transform = pFM->dom.Periodize(neighIndex,index);
        assert(transform.IsValid()); (void)transform;
        const ScalarType neighValue = values(neighIndex);
        const ScalarType valueDiff = std::max(0., value-neighValue);
        const ScalarType w = diff.Weight(data.mult)*valueDiff;
        weightSum+=w;
        neighbors.push_back({neighIndex,w});
        DiffHelper::Elem(var,diff)+=w*valueDiff;
    };
    
    const int iMax = HFM::StencilType::GetIMax(active);
    int iNeigh = iMax*HFM::StencilType::nSingleNeigh;
    const auto & forward = data.stencil.forward[iMax];
    const auto & symmetric = data.stencil.symmetric[iMax];
    
    for(const auto & diff : forward){
        if(active[iNeigh]){
            func( 1,diff);}
        ++iNeigh;
    }
    for(const auto & diff : symmetric){
        if(active[iNeigh]){
            assert(!active[iNeigh+1]);
            func( 1,diff);}
        if(active[iNeigh+1]){
            assert(!active[iNeigh]);
            func(-1,diff);}
        iNeigh+=2;
    }
    
    
    if(weightSum==0) return var;
    const ScalarType invWeight = 1./weightSum;
    for(auto & neigh : neighbors){neigh.weight*=invWeight;}
    return DiffHelper::Times(invWeight,DiffHelper::Times(data.mult,var));
}
*/

// ---------- Reverse auto diff ------------

template<typename T> auto FirstVariation<T>::
BackwardVariation(const std::vector<std::pair<IndexType,ScalarType> > & base,
				  Array<MultType, Dimension> & result) const
-> std::vector<std::pair<IndexType,ScalarType> > {
    std::vector<std::pair<IndexType,ScalarType> > sensitiveSeeds;
    const Array<ScalarType, Dimension> & values = pFM->values;
    if(result.empty()){result.dims = values.dims;
        result.resize(values.size(),DiffHelper::NullMult());}
    if(result.dims!=values.dims || result.size()!=values.size()){
        ExceptionMacro("BackwardVariation error : inconsistent input size");}
    
    // (value,index) -> weight
    typedef std::map<std::pair<ScalarType,IndexType>,ScalarType> GeoType;
    GeoType geo;
    for(auto [ind,weight] : base){
		IndexType index = ind; // Workaround for MSVC compiler bug, which forgets to discard const qualifier
        if(!pFM->dom.PeriodizeNoBase(index).IsValid()) continue;
        geo.insert({{values(index),index},weight});}

    
    while(!geo.empty()){
        const auto it = --geo.end(); // point with largest value
		const auto [value,index] = it->first;
        const ScalarType weight = it->second;
        geo.erase(it);
        
        if(weight==0) continue;
        if(pCurrentTime!=nullptr) *pCurrentTime=value;
        
        ActiveNeighborsType neighbors;
        MultType var = ValueVariation(index, neighbors);
        // Possible improvement : value(index) and value(neighbor.index) are currently accessed twice
        if(neighbors.empty()){sensitiveSeeds.push_back({index,weight});} // At seed
        
        result(index)=DiffHelper::Times(weight,var);
        // Push the weight onto the children
        for(const auto & neigh : neighbors){
            const ScalarType neighValue = pFM->values(neigh.index);
            auto it = geo.find({neighValue,neigh.index});
            if(it==geo.end()) {geo[{neighValue,neigh.index}] = neigh.weight * weight;}
            else {it->second += neigh.weight * weight;}
        }
    }
    return sensitiveSeeds;
}

// --------- Forward auto-diff ----------

template<typename T> void FirstVariation<T>::
ForwardVariation(const Array<MultType,Dimension+1> & _multVar, Array<ScalarType, Dimension+1> & result) const {
    const DiscreteType nVar = std::max(_multVar.dims[Dimension],result.dims[Dimension]);
    const IndexType dims = pFM->values.dims;
    const DeepIndexType deepDims = DeepIndex(dims,nVar);
    const DiscreteType deepSize = deepDims.Product();
    
    Array<MultType,Dimension+1> __multVar;
    const Array<MultType, Dimension+1> & multVar = _multVar.empty() ? __multVar : _multVar;
    if(multVar.empty()){__multVar.dims = deepDims; __multVar.resize(deepSize,DiffHelper::NullMult());}
    if(result.empty()){result.dims = deepDims; result.resize(deepSize,0);}
    if(result.dims != deepDims || multVar.dims!=deepDims || result.size()!=deepSize || multVar.size()!=deepSize)
        ExceptionMacro("Forward variation error : inconsistent input sizes");
    
    // Sort all values increasingly
    std::vector<std::pair<ScalarType, DiscreteType> > valIndex;
    const auto & values = pFM->values;
    valIndex.reserve(values.size());
    for(DiscreteType i=0; i<values.size(); ++i){
        valIndex.push_back({values[i],i});}
    std::sort(valIndex.begin(), valIndex.end());
    
    for(const auto & vI : valIndex){
        const IndexType index = values.Convert(vI.second);
        if(pCurrentTime!=nullptr) *pCurrentTime=vI.first;
        
        ActiveNeighborsType neighbors;
        MultType var = ValueVariation(index, neighbors);
        
        for(int k=0; k<nVar; ++k){
            const DeepIndexType dIndex = DeepIndex(index,k);
            ScalarType & res = result(dIndex);
            res+=DiffHelper::Scal(var,multVar(dIndex));
            
            for(const NeighborType & neigh : neighbors){
                const DeepIndexType dNeigh = DeepIndex(neigh.index,k);
                res+=neigh.weight*result(dNeigh);
            }
        }
    }
}


#endif /* FirstVariation_h */
