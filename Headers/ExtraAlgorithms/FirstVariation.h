// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef FirstVariation_h
#define FirstVariation_h

// Forward and backward differentiation of the value with respect to variations of the cost function
// (which is the inverse of the speed function), and of the seeds values.
// TODO : Make this work if no multiplier.
// TODO : make this exact as well in the case of second order differences/time varying speed field.

template<typename T> struct FirstVariation :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    Redeclare7Types(FromHFM,IndexType,ScalarType,IndexCRef,HFMI,DifferenceType,ActiveNeighFlagType,DiscreteType)
    Redeclare4Types(FromHFM,Traits,RecomputeType,DiscreteFlowType,PointType)
    Redeclare1Constant(FromHFM,Dimension)
    template<typename E, size_t n> using Array = typename HFM::template Array<E,n>;
    typedef std::conditional<HFM::hasMultiplier, typename HFM::MultiplierType, ScalarType> MultType;
    
    ScalarType * pCurrentTime=nullptr;
    struct NeighborType {IndexType index; ScalarType weight;};
    typedef CappedVector<NeighborType, HFM::nActiveNeigh> ActiveNeighborsType;
    MultiplierType ValueVariation(IndexCRef, ActiveNeighborsType &) const; // Elementary differentiation
    std::vector<std::pair<IndexType,ScalarType> >
    BackwardVariation(const std::vector<std::pair<IndexType,ScalarType> > &, Array<MultiplierType, Dimension> &) const;
    void ForwardVariation(const Array<MultiplierType, Dimension+1> &, Array<ScalarType, Dimension+1> &) const;
    
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*_pFM) override {pFM=_pFM; return true;}
protected:
    const HFM * pFM;
    template<size_t VMultSize=DifferenceType::multSize, typename Dummy=void> struct _DiffHelper;
    typedef _DiffHelper<> DiffHelper;
    typedef typename Array<ScalarType, Dimension+1>::IndexType DeepIndexType;
    DeepIndexType DeepIndex(IndexType index, DiscreteType k) const {
        DeepIndexType result;
        for(int i=0; i<Dimension; ++i) {result[i]=index[i];}
        result[Dimension]=k;
        return result;
    };
    template<bool b=HFM::hasMultiplier, typename Dummy=void> struct MultArrayIO;
};

// --------- Export -------

template<typename T> template<typename Dummy>
struct FirstVariation<T>::MultArrayIO<true,Dummy> {
    typedef typename FirstVariation<T>::HFM HFM;
    Redeclare2Types(FromHFM,HFMI,MultiplierType)
    Redeclare1Constant(FromHFM,Dimension)
    template<typename E,size_t n> using Array = typename T::template Array<E,n>; 
    static void Set(HFMI*that,std::string name,const Array<MultiplierType,Dimension> & a){
        that->io.SetArray(name,a);}
    static Array<MultiplierType,Dimension+1> Get(HFMI*that,std::string name){
        return that->io.template GetArray<MultiplierType,Dimension+1>(name);}
};

template<typename T> template<typename Dummy>
struct FirstVariation<T>::MultArrayIO<false,Dummy> {
    typedef typename FirstVariation<T>::HFM HFM;
    Redeclare2Types(FromHFM,HFMI,MultiplierType)
    Redeclare1Constant(FromHFM,Dimension)
    template<typename E,size_t n> using Array = typename T::template Array<E,n>; 
    static void Set(HFMI*that,std::string name,const Array<MultiplierType,Dimension> & a){}
    static Array<MultiplierType,Dimension+1> Get(HFMI*that,std::string name){return Array<MultiplierType,Dimension+1>();}
};

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
        Array<MultiplierType,Dimension> sensitivity;
        for(int i=0; i<lengths.size(); ++i){
            for(int k=0; k<lengths[i]; ++k, ++indIt, ++wIt){
                PointType p = that->stencil.Param().ADim(*indIt);
                if(that->pFM->dom.Periodize(p)[Dimension]) {
                    ExceptionMacro("Error : inspectSensitivity data points are out of range");}
                iw.push_back({that->pFM->dom.IndexFromPoint(p),*wIt});
            }
            sensitivity.clear();
            auto sensitiveSeedIndices = BackwardVariation(iw,sensitivity);
            MultArrayIO<>::Set(that,"costSensitivity_"+std::to_string(i),sensitivity);
            
            std::vector<std::pair<PointType,ScalarType> > sensitiveSeeds;
            sensitiveSeeds.reserve(sensitiveSeedIndices.size());
            for(const auto & iS : sensitiveSeedIndices){
                sensitiveSeeds.push_back({that->stencil.Param().ReDim(that->pFM->dom.PointFromIndex(iS.first)),iS.second});}
            io.SetVector("seedSensitivity_"+std::to_string(i),sensitiveSeeds);
            
            iw.clear();
        }
    }
    if(io.HasField("costVariation") || io.HasField("seedValueVariation")){
        Array<MultiplierType,Dimension+1> fwdVar;
        if(io.HasField("costVariation")) fwdVar = MultArrayIO<>::Get(that,"costVariation");
        Array<ScalarType,Dimension+1> valueVariations;
        if(io.HasField("seedValueVariation")){
            const auto seeds = io.template GetVector<PointType>("seeds");
            const auto seedVariations = io.template GetArray<ScalarType,2>("seedValueVariation");
            
            const DiscreteType nVar = seedVariations.dims[0];
            valueVariations.dims=DeepIndex(that->stencil.dims, nVar);
            valueVariations.resize(valueVariations.dims.ProductOfCoordinates(),0);
            
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
    Redeclare3Types(FromHFM,ScalarType,DifferenceType,MultiplierType);
    static MultiplierType NullMultiplier(){return MultiplierType::Constant(0);}
    static ScalarType & Elem(MultiplierType & m, const DifferenceType & diff){return m[diff.multIndex];}
    //    static void Normalize(MultiplierType & m, ScalarType r, const MultiplierType & m0){
    //        for(int i=0; i<VMultSize; ++i) m[i]*=r*m0[i];}
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
    Redeclare3Types(FromHFM,ScalarType,DifferenceType,MultiplierType);
    static MultiplierType NullMultiplier(){return 0;}
    static ScalarType & Elem(MultiplierType & m, const DifferenceType & diff){return m;}
    static MultiplierType Times(ScalarType a, MultiplierType m){return a*m;}
    //    static MultiplierType Times(MultiplierType m0, MultiplierType m1){return m0*m1;} // Same signature as above.
    static ScalarType Scal(MultiplierType m0, MultiplierType m1){return m0*m1;}
};

template<typename TTraits> template<typename Dummy>
struct FirstVariation<TTraits>::_DiffHelper<0,Dummy> {
    typedef HamiltonFastMarching<TTraits> HFM;
    Redeclare3Types(FromHFM,ScalarType,DifferenceType,MultiplierType);
    static MultiplierType NullMultiplier(){return MultiplierType();}
    static ScalarType & Elem(MultiplierType & m, const DifferenceType & diff){
        static ScalarType dummyScalar=0; return dummyScalar;}
    static MultiplierType Times(ScalarType, MultiplierType){return MultiplierType();}
    static MultiplierType Times(MultiplierType, MultiplierType){return MultiplierType();}
    static ScalarType Scal(const MultiplierType & m0, const MultiplierType & m1){return 0;}
};

// ----------- Differentiation at a single point -----------

template<typename T> auto FirstVariation<T>::
ValueVariation(IndexCRef index, ActiveNeighborsType & neighbors) const ->MultiplierType {
    assert(neighbors.empty());
    const auto & values = pFM->values;
    const ScalarType value = values(index);
    
    MultiplierType var = DiffHelper::NullMultiplier();
    ScalarType weightSum=0;
    const auto & data = pFM->pStencilData->RecomputeData(index);
    const ActiveNeighFlagType active = pFM->activeNeighs(index);
    
    auto func = [&,this](DiscreteType s, const DifferenceType & diff){
        IndexType neighIndex=index;
        for(int i=0; i<Dimension; ++i) neighIndex[i]+=s*diff.offset[i];
        const auto reversed = pFM->dom.Periodize(neighIndex);
        assert(!reversed[Dimension]);
        const ScalarType neighValue = values(neighIndex);
        const ScalarType valueDiff = std::max(0., value-neighValue);
        const ScalarType w = diff.Weight(data.mult)*valueDiff;
        weightSum+=w;
        neighbors.push_back({neighIndex,w});
        DiffHelper::Elem(var,diff)+=w*valueDiff;
    };
    
    int iNeigh=0;
    for(const auto & diff : data.stencil.forward){
        if(active[iNeigh]){
            func( 1,diff);}
        ++iNeigh;}
    for(const auto & diff : data.stencil.symmetric){
        if(active[iNeigh]){
            assert(!active[iNeigh+1]);
            func( 1,diff);}
        if(active[iNeigh+1]){
            assert(!active[iNeigh]);
            func(-1,diff);}
        iNeigh+=2;
    }
    
    int iMax=0;
    for(int i=0; i<HFM::nMaxBits; ++i){
        if(active[HFM::nNeigh+i]) iMax |= (1<<i);}
    iNeigh+=iMax*(T::nMaxForward+2*T::nMaxSymmetric);
    const auto & maxForward = data.stencil.maxForward[iMax];
    const auto & maxSymmetric = data.stencil.maxSymmetric[iMax];
    
    for(const auto & diff : maxForward){
        if(active[iNeigh]){
            func( 1,diff);}
        ++iNeigh;
    }
    for(const auto & diff : maxSymmetric){
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

// ---------- Reverse auto diff ------------

template<typename T> auto FirstVariation<T>::
BackwardVariation(const std::vector<std::pair<IndexType,ScalarType> > & base, Array<MultiplierType, Dimension> & result) const
-> std::vector<std::pair<IndexType,ScalarType> > {
    std::vector<std::pair<IndexType,ScalarType> > sensitiveSeeds;
    const Array<ScalarType, Dimension> & values = pFM->values;
    if(result.empty()){result.dims = values.dims;
        result.resize(values.size(),DiffHelper::NullMultiplier());}
    if(result.dims!=values.dims || result.size()!=values.size()){
        ExceptionMacro("BackwardVariation error : inconsistent input size");}
    
    // (value,index) -> weight
    typedef std::map<std::pair<ScalarType,IndexType>,ScalarType> GeoType;
    GeoType geo;
    for(const auto & indexWeight : base){
        IndexType index = indexWeight.first;
        if(pFM->dom.Periodize(index)[Dimension]) continue;
        geo.insert({{values(index),index},indexWeight.second});}

    
    while(!geo.empty()){
        const auto it = --geo.end(); // point with largest value
        const ScalarType value = it->first.first;
        const IndexType index = it->first.second;
        const ScalarType weight = it->second;
        geo.erase(it);
        
        if(weight==0) continue;
        if(pCurrentTime!=nullptr) *pCurrentTime=value;
        
        ActiveNeighborsType neighbors;
        MultiplierType var = ValueVariation(index, neighbors);
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
ForwardVariation(const Array<MultiplierType,Dimension+1> & _multVar, Array<ScalarType, Dimension+1> & result) const {
    const DiscreteType nVar = std::max(_multVar.dims[Dimension],result.dims[Dimension]);
    const IndexType dims = pFM->values.dims;
    const DeepIndexType deepDims = DeepIndex(dims,nVar);
    const DiscreteType deepSize = deepDims.ProductOfCoordinates();
    
    Array<MultiplierType,Dimension+1> __multVar;
    const Array<MultiplierType, Dimension+1> & multVar = _multVar.empty() ? __multVar : _multVar;
    if(multVar.empty()){__multVar.dims = deepDims; __multVar.resize(deepSize,DiffHelper::NullMultiplier());}
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
        MultiplierType var = ValueVariation(index, neighbors);
        
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
