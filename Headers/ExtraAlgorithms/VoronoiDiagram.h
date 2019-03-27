// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef VoronoiDiagram_h
#define VoronoiDiagram_h


template<typename T> struct VoronoiDiagram :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare15Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,
					 DiscreteType,ShortType,OffsetCRef,VectorType,DiscreteFlowType,
					 RecomputeType,DiscreteFlowElement,Decision,IndexDiff)
    Redeclare1Constant(HFM,Dimension)
    
    enum StoppingCriterionEnum {kVoronoiNone,kVoronoiRegionsMeeting,kVoronoiOppositesMeeting};
    StoppingCriterionEnum stoppingCriterion = kVoronoiNone;
    IndexType stoppingIndex0=IndexType::Constant(-1), stoppingIndex1=IndexType::Constant(-1); // Set when voronoi regions meet
	typedef ShortType VoronoiFlagType;
    typename HFM::template Array<VoronoiFlagType,Dimension> voronoiFlags;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) override;
    const HFM * pFM=nullptr;
    template<bool hasBundle_,typename Dummy> struct OppositeIndexFunctor_;
	typedef OppositeIndexFunctor_<HFM::hasBundle, void> OppositeIndexFunctor;
};

template<typename T> void VoronoiDiagram<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(!(io.HasField("seedFlags") || (HFM::hasBundle && io.HasField("seedFlags_Unoriented")))) return;
    voronoiFlags.dims = that->stencil.dims;
    voronoiFlags.resize(voronoiFlags.dims.Product(),-1);
    
    if(io.HasField("seedFlags")){
        const auto seedFlags = io.template GetVector<ScalarType>("seedFlags");
        const auto seeds = io.template GetVector<PointType>("seeds");
        if(seedFlags.size()!=seeds.size()){
            ExceptionMacro("Error : Number of seeds (" <<seeds.size()<<
                           ") is distinct from the number of seedFlags(" <<seedFlags.size() << ").\n");}
        for(size_t i=0; i<seeds.size(); ++i){
            IndexType index = that->pFM->dom.IndexFromPoint(that->stencil.Param().ADim(seeds[i]));
            if(!that->pFM->dom.Periodize(index,index).IsValid()) ExceptionMacro("Error : seed " << seeds[i] << " is out of range.\n");
            voronoiFlags(index) = (VoronoiFlagType)seedFlags[i];

        }
    }
    
    if(HFM::hasBundle && io.HasField("seedFlags_Unoriented")) {
        typedef typename HFMI::SpecializationsDefault HFMIS;
        const auto uSeedFlags = io.template GetVector<ScalarType>("seedFlags_Unoriented");
        const auto uSeeds = io.template GetVector<typename HFMIS::UnorientedPointType>("seeds_Unoriented");
        if(uSeeds.size()!=uSeedFlags.size()){
            ExceptionMacro("Error : Number of unoriented seeds (" <<uSeeds.size()<<
                           ") is distinct from the number of unoriented seedFlags(" <<uSeedFlags.size() << ").\n");}
        std::vector<PointType> equiv;
        for(int i=0; i<uSeeds.size(); ++i){
            HFMIS::PadAdimEquiv(that,uSeeds[i],equiv);
            for(const PointType & p : equiv){
                IndexType index = that->pFM->dom.IndexFromPoint(p);
                if(!that->pFM->dom.Periodize(index,index).IsValid()) ExceptionMacro("Error : unoriented seed, of index " << i << ", is out of range.\n");
                voronoiFlags(index) = (VoronoiFlagType)uSeedFlags[i];
            }
        }
    }
    
    const std::string s = io.GetString("voronoiStoppingCriterion","None");
    if(s=="RegionsMeeting") stoppingCriterion=StoppingCriterionEnum::kVoronoiRegionsMeeting;
    else if(HFM::hasBundle && s=="OppositesMeeting") stoppingCriterion=StoppingCriterionEnum::kVoronoiOppositesMeeting;
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
struct VoronoiDiagram<Traits>::OppositeIndexFunctor_<true,Dummy>{
    Redeclare2Types(Traits,IndexType,DiscreteType);
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
struct VoronoiDiagram<Traits>::OppositeIndexFunctor_<false,Dummy>{
    Redeclare1Type(Traits,IndexType);
    IndexType operator()(IndexType,IndexType) const {assert(false); return IndexType();}
};

template<typename T> int VoronoiDiagram<T>::
PostProcessWithRecompute(IndexCRef index, const RecomputeType &, const DiscreteFlowType & flow) {
    
    VoronoiFlagType & indexFlag = voronoiFlags(index);
    ScalarType maxWeight=0;
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight<=maxWeight) continue;
        maxWeight=fl.weight;
        IndexType neigh = index + IndexDiff::CastCoordinates(fl.offset);
        const auto transform = pFM->dom.Periodize(neigh,index);
        assert(transform.IsValid()); (void)transform;
        indexFlag = voronoiFlags(neigh);
    }
    
    if(stoppingCriterion==StoppingCriterionEnum::kVoronoiRegionsMeeting){ // Note: this would not work for a Dijkstra
        for(DiscreteType i=0; i<Dimension; ++i){ // Look for a neighbor with a distinct index
            for(DiscreteType eps=-1; eps<=1; ++eps){
                IndexType neigh=index;
                neigh[i]+=eps;
                const auto transform = pFM->dom.Periodize(neigh,index);
                if(!transform.IsValid()) continue;
                const VoronoiFlagType flag = voronoiFlags(neigh);
                if(flag==-1 || flag==indexFlag) continue;
                stoppingIndex0=index;
                stoppingIndex1=neigh;
                return Decision::kTerminate;
            }
        }
        /*
        ShortType candidateFlag=-1;
        IndexType candidateIndex;
        for(const DiscreteFlowElement & fl : flow){
            if(fl.weight==0) continue;
            IndexType neigh = index; for(int i=0; i<Dimension; ++i) neigh[i]+=fl.offset[i];
            const auto reversed = pFM->dom.Periodize(neigh);
            assert(!reversed[Dimension]); (void)reversed;
            ShortType flag = voronoiFlags(neigh);
            if(candidateFlag==-1) {candidateFlag=flag;candidateIndex=neigh;}
            if(flag!=candidateFlag){
                stoppingIndex0=candidateIndex;
                stoppingIndex1=neigh;
                return Decision::kTerminate;
            }
        }*/
    }
    
    if(stoppingCriterion==StoppingCriterionEnum::kVoronoiOppositesMeeting){ // Note : only applies to R^n x S^d domains, d in {1,2}
        const IndexType oppositeIndex = OppositeIndexFunctor()(index,pFM->stencilData.dims);
        const VoronoiFlagType oppositeFlag = voronoiFlags(oppositeIndex);
        if(indexFlag!=-1 && oppositeFlag!=-1 && indexFlag!=oppositeFlag){
            stoppingIndex0=index; stoppingIndex1=oppositeIndex;
            return Decision::kTerminate;
        }
    }
    
    return 0;
}

#endif /* VoronoiDiagram_h */
