// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef EuclideanPathLength_h
#define EuclideanPathLength_h

/*
 // TODO : (why not) geodesic from stopping point.
 */

template<typename T> struct EuclideanPathLength :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare13Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,
					 DiscreteType,OffsetCRef,VectorType,DiscreteFlowType,RecomputeType,
					 DiscreteFlowElement,Decision)
    Redeclare1Constant(HFM,Dimension)
    
    ScalarType stopAtEuclideanLength = Traits::Infinity();
    IndexType stoppingIndex = IndexType::Constant(-1); // Output
    std::unique_ptr<typename HFM::template DataSource<VectorType> > euclideanScale;
    typename HFM::template Array<ScalarType,Dimension> euclideanLengths;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) override;
    const HFM * pFM=nullptr;
};

template<typename T> void EuclideanPathLength<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(io.HasField("euclideanScale")) {
        euclideanScale = that->template GetField<VectorType>("euclideanScale");
        if(io.HasField("stopAtEuclideanLength")){
            stopAtEuclideanLength = io.template Get<ScalarType>("stopAtEuclideanLength");}
    }
}

template<typename T> void EuclideanPathLength<T>::
Finally(HFMI*that){
    auto & io = that->io;
    if(io.template Get<ScalarType>("exportEuclideanLengths",1)) {
        io.SetArray("euclideanLengths",euclideanLengths);}
    if(stoppingIndex!=IndexType::Constant(-1)){
        io.template Set<PointType>("euclideanLength_stoppingPoint",that->stencil.Param().ReDim(pFM->dom.PointFromIndex(stoppingIndex)));
        if(io.template Get<ScalarType>("euclideanLength_exportGeodesicFromStoppingPoint",1.)){
            that->ExportGeodesics("euclideanLength",{that->pFM->dom.PointFromIndex(stoppingIndex)});}
    }
}

template<typename T> bool EuclideanPathLength<T>::
ImplementIn(HFM*_pFM){
    if(euclideanScale==nullptr) return false;
    pFM=_pFM;
    _pFM->extras.postProcessWithRecompute.push_back(this);
    euclideanLengths.dims=pFM->values.dims;
    euclideanLengths.resize(euclideanLengths.dims.Product(),-1);
    return true;
}

template<typename T> int EuclideanPathLength<T>::
PostProcessWithRecompute(IndexCRef index, const RecomputeType &, const DiscreteFlowType & flow){
    assert(pFM!=nullptr); assert(!euclideanLengths.empty());
    VectorType flowSum=VectorType::Constant(0);
    ScalarType lengthSum=0, weightSum=0;
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight==0) continue;
        weightSum+=fl.weight;
        flowSum+=fl.weight * VectorType::CastCoordinates(fl.offset);
        IndexType neigh = index; for(int i=0; i<Dimension; ++i) neigh[i]+=fl.offset[i];
        const auto transform = pFM->dom.Periodize(neigh,index);
        assert(transform.IsValid()); (void)transform;
        lengthSum+=fl.weight * euclideanLengths(neigh);
    }
    if(weightSum>0){
        flowSum/=weightSum;
        lengthSum/=weightSum;
        const VectorType scale = (*euclideanScale)(index);
        for(int i=0; i<Dimension; ++i) flowSum[i]*=scale[i];
        const ScalarType length = flowSum.Norm() + lengthSum;
        euclideanLengths(index) = length;
        if(length>stopAtEuclideanLength) {
            stoppingIndex = index;
            return Decision::kTerminate;}
    } else {euclideanLengths(index) = 0;}
    return 0;
}

#endif /* EuclideanPathLength_h */
