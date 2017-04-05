// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HFMInterface_hxx
#define HFMInterface_hxx

// ---- Default stencil setup ------

template<typename T> template<typename Dummy> struct
HFMInterface<T>::SpecializationsDefault<true,Dummy> {
    typedef HFMInterface<T> HFMI;
    typedef typename HFMI::HFM HFM;
    Redeclare6Types(FromHFM,ScalarType,Traits,MultiplierType,PointType,IndexType,DiscreteType)
    Redeclare2Constants(FromTraits,mathPi,Dimension)
    
    static const int DimDep =Traits::nStencilDependencies, DimIndep = Dimension-DimDep;
    template<typename E> using DataSource = typename HFM::template DataSource<E>;
    template<typename E> using DataSource_Value = typename HFMI::template DataSource_Value<E>;
    template<typename E> using DataSource_Array = typename HFMI::template DataSource_Array<E>;
    template<typename E> using DataSource_Indep = typename HFMI::template DataSource_Indep<E>;
    template<typename E> using DataSource_Dep = typename HFMI::template DataSource_Dep<E>;
    template<typename E> static std::unique_ptr<DataSource<E> > GetField(std::string name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<E>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>(io.template Get<E>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<E, Dimension>(name)));}
        else if(nDims==DimIndep) {  return ResultType(new DataSource_Indep<E>(io.template GetArray<E, DimIndep>(name)));}
        else if(nDims==DimDep) {    return ResultType(new DataSource_Dep<E>(io.template GetArray<E, DimDep>(name)));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.\n");}
    }
    template<typename E> static std::unique_ptr<DataSource<E> > GetIntegralField(std::string name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<ScalarType>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>((E)io.template Get<ScalarType>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<ScalarType, Dimension>(name).template Cast<E>()));}
        else if(nDims==DimIndep) {  return ResultType(new DataSource_Indep<E>(io.template GetArray<ScalarType, DimIndep>(name).template Cast<E>()));}
        else if(nDims==DimDep) {    return ResultType(new DataSource_Dep<E>(io.template GetArray<ScalarType, DimDep>(name).template Cast<E>()));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.");}
    }
    typedef std::array<ScalarType,DimIndep==0 ? 1 : DimIndep> UnorientedPointType;
    static PointType PadUnoriented(const UnorientedPointType & p){
        PointType q=PointType::Constant(0);
        auto depIt = Traits::stencilDependencies.begin();
        auto pIt = p.begin();
        for(int i=0; i<Dimension; ++i) {
            if(depIt!=Traits::stencilDependencies.end() && i==*depIt) ++depIt;
            else {q[i]=*pIt; ++pIt;}
        }
        assert(DimIndep==0 || pIt==p.end());
        return q;
    }
    static void EquivalentPoints(HFMI*that,PointType p, std::vector<PointType> & result){
        const IndexType dims = that->stencil.dims;
        const auto & dom = that->pFM->dom;
        IndexType index=IndexType::Constant(0);
        const auto dep = Traits::stencilDependencies;
        while(true){
            const PointType r = dom.PointFromIndex(index);
            for(DiscreteType i : dep){p[i]=r[i];}
            result.push_back(p);

            // Next index
            int k;
            for(k=0; k<Traits::nStencilDependencies; ++k){
                if(index[dep[k]]<dims[dep[k]]-1) {++index[dep[k]]; break;}
                else {index[dep[k]]=0;}
            }
            if(k==Traits::nStencilDependencies) break;
        }
    }
};

template<typename T> template<typename Dummy> struct
HFMInterface<T>::SpecializationsDefault<false,Dummy> {
    typedef HFMInterface<T> HFMI;
    typedef typename HFMI::HFM HFM;
    Redeclare2Types(FromHFM,ScalarType,PointType)
    Redeclare1Constant(FromHFM,Dimension)
    
    template<typename E> using DataSource = typename HFM::template DataSource<E>;
    template<typename E> using DataSource_Value = typename HFMI::template DataSource_Value<E>;
    template<typename E> using DataSource_Array = typename HFMI::template DataSource_Array<E>;
    template<typename E> static std::unique_ptr<DataSource<E> > GetField(std::string name, HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<E>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>(io.template Get<E>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<E, Dimension>(name)));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.\n");}
    }
    template<typename E> static std::unique_ptr<DataSource<E> > GetIntegralField(std::string name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<ScalarType>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>((E)io.template Get<ScalarType>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<ScalarType, Dimension>(name).template Cast<E>()));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.");}
    }
    static void EquivalentPoints(HFMI*,PointType, std::vector<PointType> &){assert(false);};
    typedef PointType UnorientedPointType;
    static UnorientedPointType PadUnoriented(const PointType & p){return p;}
};


// ------- Data sources ------
template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Value : DataSource<E> {
    typedef typename DataSource<E>::ReturnType ReturnType;
    E value;
    DataSource_Value(E _value):value(std::move(_value)){};
    virtual bool CheckDims(const IndexType &) const {return true;}
    virtual ReturnType operator()(const IndexType &) const {return value;}
};

template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Array : DataSource<E> {
    typedef typename DataSource<E>::ReturnType ReturnType;
    Array<E,Dimension> values;
    DataSource_Array(Array<E, Dimension> && _values):values(std::move(_values)){};
    virtual bool CheckDims(const IndexType & dims) const {return dims==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(index);}
};


template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Indep<E,true> : DataSource<E> {
    typedef typename HFMInterface<T>::HFM HFM;
    typedef typename HFM::template DataSource<E>::ReturnType ReturnType;
    Redeclare2Types(FromHFM,Traits,IndexType)
    Redeclare1Constant(FromHFM,Dimension)
    typedef HFMInterface<T>::Array<E, Dimension-Traits::nStencilDependencies> ArrayType;
    ArrayType values;
    DataSource_Indep(ArrayType && _values):values(std::move(_values)){};
    typedef typename ArrayType::IndexType ShortIndexType;
    ShortIndexType Convert(const IndexType & index) const {
        ShortIndexType result;
        int i=0;
        for(int j=0; j<Dimension; ++j){
            if(i<Traits::nStencilDependencies && Traits::stencilDependencies[i]==j){++i;}
            else result[j-i]=index[j];
        }
        return result;
    }
    virtual bool CheckDims(const IndexType & dims) const {return Convert(dims)==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(Convert(index));}
};

template<typename T> template<typename E> struct // Dummy specialization
HFMInterface<T>::DataSource_Indep<E,false> {};

template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Dep<E,true> : DataSource<E> {
    typedef typename HFMInterface<T>::HFM HFM;
    typedef typename HFM::template DataSource<E>::ReturnType ReturnType;
    Redeclare2Types(FromHFM,Traits,IndexType)
    Redeclare1Constant(FromHFM,Dimension)
    typedef HFMInterface<T>::Array<E, Traits::nStencilDependencies> ArrayType;
    ArrayType values;
    DataSource_Dep(ArrayType && _values):values(std::move(_values)){};
    typedef typename ArrayType::IndexType ShortIndexType;
    ShortIndexType Convert(const IndexType & index) const {
        ShortIndexType result;
        for(int i=0; i<Traits::nStencilDependencies; ++i){
            result[i]=index[Traits::stencilDependencies[i]];}
        return result;
    }
    virtual bool CheckDims(const IndexType & dims) const {return Convert(dims)==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(Convert(index));}
};

template<typename T> template<typename E> struct // Dummy specialization
HFMInterface<T>::DataSource_Dep<E,false> {};

template<typename T> template<typename E>
struct HFMInterface<T>::TimeDependentSource : DataSource<E> {
    typedef T Traits;
    const ScalarType * pCurrentTime=nullptr;
    typedef std::pair<ScalarType,std::unique_ptr<DataSource<E> > > InterpolationPoint;
    std::vector<InterpolationPoint> interpolationData;
    typedef typename DataSource<E>::ReturnType ReturnType;
    virtual bool CheckDims(const IndexType & index) const {
        // Also checks that times are sorted and pointers well assigned
        if(interpolationData.empty()) return false;
        ScalarType lastTime = -T::Infinity();
        for(const auto & timeSource : interpolationData){
            if(timeSource.first<=lastTime) return false;
            lastTime=timeSource.first;
            if(timeSource.second==nullptr) return false;
            if(!timeSource.second->CheckDims(index)) return false;
        }
        return true;
    };
    virtual ReturnType operator()(const IndexType & index) const {
        assert(pCurrentTime!=nullptr);
        const ScalarType t=*pCurrentTime; 
        if(t<=interpolationData.front().first) return (*interpolationData.front().second)(index);
        if(t>=interpolationData.back().first) return (*interpolationData.back().second)(index);
        const auto it1 = std::lower_bound(interpolationData.begin(), interpolationData.end(), *pCurrentTime,
                                          [](const InterpolationPoint & a, ScalarType b)->bool{return a.first<b;});
        const auto it0 = it1-1;
        const ScalarType t0 = it0->first, t1=it1->first;
        const E r0 = (*it0->second)(index), r1=(*it1->second)(index);
        return r0 + ((t-t0)/(t1-t0))*(r1-r0);
    };
};

// ---- Getting fields -------

template<typename T> template<typename E> auto
HFMInterface<T>::GetField(std::string s, bool mayDependOnTime) -> std::unique_ptr<DataSource<E> > {
    if(io.HasField(s)){
    std::unique_ptr<DataSource<E> > result = SpecializationsDefault<>::template GetField<E>(s,this);
    if(!result->CheckDims(stencil.dims)){
        ExceptionMacro("Error : Field " << s << " has inconsistent dimensions");}
        return std::move(result);
    } else if(mayDependOnTime && io.HasField(s+"_times")) {
        typedef TimeDependentSource<E> ResultType;
        std::unique_ptr<ResultType> pResult(new ResultType);
        pResult->pCurrentTime = &(pTime->currentTime);
        pTime->currentTime = -Traits::Infinity();
        std::vector<ScalarType> times = io.template GetVector<ScalarType>(s+"_times");
        for(int i=1; i<times.size(); ++i){
            if(times[i-1]>=times[i]) ExceptionMacro("Time dependent field error : times are not strictly increasing.");}
        pResult->interpolationData.reserve(times.size());
        for(int i=0; i<times.size(); ++i){
            pResult->interpolationData.push_back({times[i], GetField<E>(s+"_"+std::to_string(i))});}
        return std::move(pResult);
    } else ExceptionMacro("Error : Field " << s << " not found.\n");
    
}

template<typename T> template<typename E> auto
HFMInterface<T>::GetIntegralField(std::string s) -> std::unique_ptr<DataSource<E> > {
    std::unique_ptr<DataSource<E> > result = SpecializationsDefault<>::template GetIntegralField<E>(s,this);
    if(!result->CheckDims(stencil.dims)){
        ExceptionMacro("Error : Field " << s << " has inconsistent dimensions.\n");}
    return std::move(result);}


// ---------- Running  ---------
template<typename T> void HFMInterface<T>::
Run() {
    Run_SetupIO();
    Run_SetupStencil();
    Run_SetupSolver();
    Run_SetupExtraAlgorithms();
    if(Run_RunSolver()) return;
    Run_ExtractGeodesics();
    Run_ExportData();
}

template<typename T> void HFMInterface<T>::
Run_SetupIO() {
    // Setup input array ordering policy
    io.verbosity = io.Get<ScalarType>("verbosity",io.verbosity);
    io.arrayOrdering = enumFromString<ArrayOrdering>(io.GetString("arrayOrdering",enumToRealString(io.arrayOrdering)));
    pTime = std::unique_ptr<TimeDependentFields<T> >(new TimeDependentFields<T>);
    pTime->Setup(this);
}

template<typename T> void HFMInterface<T>::
Run_SetupStencil() {
    // Setup stencils (hence metric)
    stencil.Setup(this);
    
    const clock_t top = clock();
    pFM=std::unique_ptr<HFM>(new HFM(std::move(_pStencil)));
    io.Set<ScalarType>("StencilCPUTime",ScalarType(clock()-top)/CLOCKS_PER_SEC);
}

template<typename T> void HFMInterface<T>::
Run_SetupSolver() {
    // ------- Some exports that are independent of the fast marching results -------
    io.Set<ScalarType>("MaxStencilWidth",pFM->MaxStencilWidth());
    pFM->sndOrder = (bool)io.template Get<ScalarType>("sndOrder",0.);

    if(io.HasField("getStencils")) {
        const auto & pts = io.GetVector<PointType>("getStencils");
        std::ostringstream oss;
        oss << "{";
        for(const PointType & p : pts){
            IndexType index = pFM->dom.IndexFromPoint(stencil.Param().ADim(p));
            if(pFM->dom.Periodize(index)[Dimension]){
                WarnMsg() << "GetStencil error : point " << p << "is out of range.\n";
                continue;}
            typename HFM::StencilType indStencil;
            stencil.SetStencil(index,indStencil);
            oss << "{" << index << "," << indStencil << "},";
        }
        oss << "}";
        io.SetString("stencils", oss.str());
    }
    
    if(io.HasField("pointToIndex")){
        auto pts = io.GetVector<PointType>("pointToIndex");
        for(PointType & p : pts) p=PointType::CastCoordinates(pFM->dom.IndexFromPoint(stencil.Param().ADim(p)));
        io.SetVector("indexFromPoint", pts);
    }
    if(io.HasField("indexToPoint")){
        auto pts = io.GetVector<PointType>("indexToPoint");
        for(PointType & p : pts) p=stencil.Param().ReDim(pFM->dom.PointFromIndex(IndexType::CastCoordinates(p)));
        io.SetVector("pointFromIndex", pts);
    }

    
    // Get the seeds
    {
        std::vector<PointType> seedPoints;
        std::vector<ScalarType> seedValues;
        if(io.HasField("seeds")){
            seedPoints = io.GetVector<PointType>("seeds");
            for(auto & p : seedPoints) p=stencil.Param().ADim(p);
            if(io.HasField("seedValues")) seedValues = io.GetVector<ScalarType>("seedValues");
            else seedValues.resize(seedPoints.size(),0.);
        }
        
        if(HFM::hasMultiplier && io.HasField("seeds_Unoriented")){
            typedef typename SpecializationsDefault<>::UnorientedPointType UnorientedPointType;
            const auto _uPoints = io.GetVector<UnorientedPointType>("seeds_Unoriented");
            std::vector<PointType> uPoints; uPoints.reserve(_uPoints.size());
            for(const auto & p : _uPoints) uPoints.push_back(SpecializationsDefault<>::PadUnoriented(p));
                
            std::vector<ScalarType> uValues;
            if(io.HasField("seedValues_Unoriented")) uValues = io.GetVector<ScalarType>("seedValues_Unoriented");
            else uValues.resize(uPoints.size(),0.);
            if(uValues.size()!=uPoints.size()) ExceptionMacro("Error : Inconsistent size of seedValues_Unoriented.");
            
            std::vector<PointType> equiv;
            for(int i=0; i<uPoints.size(); ++i){
                equiv.clear();
                SpecializationsDefault<>::EquivalentPoints(this,stencil.Param().ADim(uPoints[i]),equiv);
                seedPoints.insert(seedPoints.end(),equiv.begin(),equiv.end());
                seedValues.resize(seedValues.size()+equiv.size(),uValues[i]);
            }
        }
        
        for(int i=0; i<seedPoints.size(); ++i){
            auto seedIndex=pFM->dom.IndexFromPoint(seedPoints[i]);
            if(pFM->dom.Periodize(seedIndex)[Dimension]){
                WarnMsg() << "Error : seed " << seedPoints[i] << " is out of range.\n";
                continue;}
            pFM->seeds.insert({seedIndex,seedValues[i]});}
        if(pFM->seeds.empty() && !seedPoints.empty())
            ExceptionMacro("Error : seeds incorrectly set");
    }
    
    // Reinserting data from previous run
    if(io.HasField("values")){
        pFM->values = io.GetArray<ScalarType, Dimension>("values");}
    
    typedef ActiveNeighFlagType ActiveFlag;
    if(io.HasField("activeNeighs")){
        pFM->activeNeighs = io.GetArray<ScalarType, Dimension>("activeNeighs")
        .template Transform<>([](ScalarType a)->ActiveFlag {
            return ActiveFlag(long(a));});}
}

template<typename T> template<typename Alg> Alg* HFMInterface<T>::
SetupSingleAlgorithm(){
    auto pAlg = std::unique_ptr<Alg>(new Alg);
    pAlg->Setup(this);
    if(! pAlg->ImplementIn(pFM.get())) return nullptr;
    Alg * result = pAlg.get();
    extras.push_back(std::move(pAlg));
    return result;
}

template<typename T> bool HFMInterface<T>::
Run_RunSolver() {
    if(!pFM->seeds.empty()){
        const clock_t top = clock();
        pFM->Run();
        const clock_t elapsed = clock()-top;
        const ScalarType FMCPUTime = ScalarType(elapsed)/CLOCKS_PER_SEC;
        io.Set<ScalarType>("FMCPUTime",FMCPUTime);
        if(io.verbosity>=1) Msg() << "Fast marching solver completed in " << FMCPUTime << " s.\n";
    } else {
        if(!io.HasField("values") || !io.HasField("activeNeighs")){
            WarnMsg() << "No seeds to run fast marching solver,"
            << " and missing values or activeNeighs to compute geodesics.\n";
            return true;
        }
    }
    return false;
}

template<typename T> void HFMInterface<T>::
ExportGeodesics(std::string suffix, const std::vector<PointType> & tips){
    if(pGeodesicSolver==nullptr){
        typedef std::unique_ptr<GeodesicSolverInterface> GeoSolverPtr;
        std::string geodesicSolverType = io.GetString("geodesicSolver", "Discrete");
        if(geodesicSolverType=="Discrete") {pGeodesicSolver = GeoSolverPtr(new GeodesicDiscreteSolver<Traits>(*pFM));}
        else if(geodesicSolverType=="ODE") {pGeodesicSolver = GeoSolverPtr(new GeodesicODESolver<Traits>(*pFM));}
        else {ExceptionMacro("Error : Unrecognized geodesic solver " << geodesicSolverType << ".\n");}
        pGeodesicSolver->Setup(this);
    }
    
    const auto & geodesics = pGeodesicSolver->Run(this,tips);

    std::vector<ScalarType> geodesicLengths;
    for(const auto & geo : geodesics) {geodesicLengths.push_back(geo.size());}
    std::vector<PointType> geodesicPoints;
    geodesicPoints.reserve((size_t)std::accumulate(geodesicLengths.begin(), geodesicLengths.end(), 0.));
    for(const auto & geo : geodesics) for(const PointType & p : geo) geodesicPoints.push_back(stencil.Param().ReDim(p));
    
    if(!suffix.empty()) suffix = "_"+suffix;
    io.SetVector("geodesicPoints"+suffix, geodesicPoints);
    io.SetVector("geodesicLengths"+suffix, geodesicLengths);
}

template<typename T> void HFMInterface<T>::
Run_ExtractGeodesics() {
    if(io.HasField("tips")){
        std::vector<PointType> tips = io.GetVector<PointType>("tips");
        for(PointType & tip : tips) tip = stencil.Param().ADim(tip);
        ExportGeodesics("",tips);
    }
    if(HFM::hasMultiplier && io.HasField("tips_Unoriented")){
        typedef typename SpecializationsDefault<>::UnorientedPointType UnorientedPointType;
        const auto _indepTips = io.GetVector<UnorientedPointType>("tips_Unoriented");
        std::vector<PointType> indepTips; indepTips.reserve(_indepTips.size());
        for(const auto & p : _indepTips) indepTips.push_back(SpecializationsDefault<>::PadUnoriented(p));
            
        std::vector<PointType> tips;
        std::vector<PointType> equiv;
        for(const PointType & indepTip : indepTips){
            PointType p = stencil.Param().ADim(indepTip);
            equiv.clear();
            SpecializationsDefault<>::EquivalentPoints(this,p,equiv);
            ScalarType valMin=Traits::Infinity();
            for(const PointType & q : equiv){
                PointType r=q;
                if(pFM->dom.Periodize(r)[Dimension]) continue; // Out of range
                ScalarType val =pFM->values(pFM->dom.IndexFromPoint(r));
                if(valMin<val) continue; // Too large value
                p=q; valMin=val;
            }
            tips.push_back(p);
        }
        ExportGeodesics("Unoriented", tips);
    }
}


template<typename T> void HFMInterface<T>::
Run_ExportData() {
    if(io.Get<ScalarType>("exportValues",0.)) {
        io.SetArray("values", pFM->values);}
    
    if(io.Get<ScalarType>("exportActiveNeighs",0.)) {
        io.SetArray("activeNeighs",pFM->activeNeighs.template Transform<>(
    [](ActiveNeighFlagType a)->ScalarType{return ScalarType(a.to_ulong());}
    ));}
    
    if(io.Get<ScalarType>("exportGeodesicFlow",0.)) {
        Array<VectorType, Dimension> flow;
        flow.dims=pFM->values.dims;
        flow.resize(pFM->values.size());
        for(DiscreteType i=0; i<flow.size(); ++i){
            flow[i] = pFM->GeodesicFlow(flow.Convert(i)).flow;} //TODO : ReDim ?
        io.SetArray<VectorType>("geodesicFlow", flow);
    }
    
    for(const auto & pAlg : extras) pAlg->Finally(this);
}

#endif /* HFMInterface_hxx */
