// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef VoronoiDiagram_h
#define VoronoiDiagram_h


template<typename T> struct VoronoiDiagram :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare8Types(FromHFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,DiscreteType,ShortType)
    Redeclare6Types(FromHFM,OffsetCRef,VectorType,DiscreteFlowType,RecomputeType,DiscreteFlowElement,Decision)
    Redeclare1Constant(FromHFM,Dimension)
    
    enum StoppingCriterionEnum {kVoronoiNone,kVoronoiRegionsMeeting,kVoronoiOppositesMeeting};
    StoppingCriterionEnum stoppingCriterion = kVoronoiNone;
    IndexType stoppingIndex0=IndexType::Constant(-1), stoppingIndex1=IndexType::Constant(-1); // Set when voronoi regions meet
    typename HFM::template Array<ShortType,Dimension> voronoiFlags;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) override;
    const HFM * pFM=nullptr;
    template<bool b=HFM::hasMultiplier,typename Dummy=void> struct OppositeIndex;
};

template<typename T> void VoronoiDiagram<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(!io.HasField("seedFlags")) return;
    const auto seedFlags = io.template GetVector<ScalarType>("seedFlags");
    const auto seeds = io.template GetVector<PointType>("seeds");
    if(seedFlags.size()!=seeds.size()){
        ExceptionMacro("Error : Number of seeds (" <<seeds.size()<<
                       ") distinct from the number of seedFlags(" <<seedFlags.size() << ").\n");}
    voronoiFlags.dims = that->stencil.dims;
    voronoiFlags.resize(voronoiFlags.dims.ProductOfCoordinates(),-1);
    for(int i=0; i<seeds.size(); ++i){
        voronoiFlags(that->pFM->dom.IndexFromPoint(that->stencil.Param().ADim(seeds[i]))) = seedFlags[i];
    }
    const std::string s = io.GetString("voronoiStoppingCriterion","None");
    if(s=="RegionsMeeting") stoppingCriterion=StoppingCriterionEnum::kVoronoiRegionsMeeting;
    else if(HFM::hasMultiplier && s=="OppositesMeeting") stoppingCriterion=StoppingCriterionEnum::kVoronoiOppositesMeeting;
    else if(s!="None") ExceptionMacro("Voronoi diagram error : unrecognized stopping criterion");
}

template<typename T> bool VoronoiDiagram<T>::
ImplementIn(HFM*_pFM){
    if(voronoiFlags.empty()) return false;
    _pFM->extras.postProcessWithRecompute.push_back(this);
    pFM=_pFM;
    assert(voronoiFlags.dims==pFM->values.dims);
    return true;
}

template<typename T> void VoronoiDiagram<T>::
Finally(HFMI*that){
    auto & io = that->io;
    if(io.template Get<ScalarType>("exportVoronoiFlags",1)) {
        io.SetArray("voronoiFlags",voronoiFlags.template Cast<ScalarType>());}
    if(stoppingIndex0!=IndexType::Constant(-1)){
        const PointType p0 = pFM->dom.PointFromIndex(stoppingIndex0), p1 = pFM->dom.PointFromIndex(stoppingIndex1);
        io.template Set<PointType>("voronoiDiagram_meetingPoint0",that->stencil.Param().ReDim(p0));
        io.template Set<PointType>("voronoiDiagram_meetingPoint1",that->stencil.Param().ReDim(p1));
        
        if(io.template Get<ScalarType>("voronoiDiagram_exportGeodesicFromMeetingPoint",1)){
            that->ExportGeodesics("voronoiDiagram",{p0,p1});}
    
    }
}

template<typename Traits> template<typename Dummy>
struct VoronoiDiagram<Traits>::OppositeIndex<true,Dummy>{
    Redeclare2Types(FromTraits,IndexType,DiscreteType);
    IndexType operator()(IndexType index,IndexType dims) const {  // Only works for R^n x S^d, d in {1,2}
            IndexType result=index;
            const auto & dep = Traits::stencilDependencies;
            const size_t nDep = dep.size();
            if(nDep>2){assert(false);}
            if(nDep>=1){DiscreteType i = dep.back(); result[i] = PosMod(result[i]+dims[i]/2, dims[i]);}
            if(nDep==2){DiscreteType i = dep.front(); result[i] = dims[i]-1-result[i];}
            return result;
    }
};

template<typename Traits> template<typename Dummy>
struct VoronoiDiagram<Traits>::OppositeIndex<false,Dummy>{
    Redeclare1Type(FromTraits,IndexType);
    IndexType operator()(IndexType,IndexType) const {assert(false); return IndexType();}
};

template<typename T> int VoronoiDiagram<T>::
PostProcessWithRecompute(IndexCRef index, const RecomputeType &, const DiscreteFlowType & flow) {
    
    ShortType & indexFlag = voronoiFlags(index);
    ScalarType maxWeight=0;
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight<=maxWeight) continue;
        maxWeight=fl.weight;
        IndexType neigh = index; for(int i=0; i<Dimension; ++i) neigh[i]+=fl.offset[i];
        const auto reversed = pFM->dom.Periodize(neigh); assert(!reversed[Dimension]);
        indexFlag = voronoiFlags(neigh);
    }
    
    if(stoppingCriterion==StoppingCriterionEnum::kVoronoiRegionsMeeting){ // Note: this would not work for a Dijkstra
        ShortType candidateFlag=-1;
        IndexType candidateIndex;
        for(const DiscreteFlowElement & fl : flow){
            if(fl.weight==0) continue;
            IndexType neigh = index; for(int i=0; i<Dimension; ++i) neigh[i]+=fl.offset[i];
            const auto reversed = pFM->dom.Periodize(neigh); assert(!reversed[Dimension]);
            ShortType flag = voronoiFlags(neigh);
            if(candidateFlag==-1) {candidateFlag=flag;candidateIndex=neigh;}
            if(flag!=candidateFlag){
                stoppingIndex0=candidateIndex;
                stoppingIndex1=neigh;
                return Decision::kTerminate;
            }
        }
    }
    if(stoppingCriterion==StoppingCriterionEnum::kVoronoiOppositesMeeting){ // Note : only applies to R^n x S^d domains, d in {1,2}
        const IndexType oppositeIndex = OppositeIndex<>()(index,pFM->pStencilData->dims);
        const ShortType oppositeFlag = voronoiFlags(oppositeIndex);
        if(indexFlag!=-1 && oppositeFlag!=-1 && indexFlag!=oppositeFlag){
            stoppingIndex0=index; stoppingIndex1=oppositeIndex;
            return Decision::kTerminate;
        }
    }
    
    return 0;
}

#endif /* VoronoiDiagram_h */
