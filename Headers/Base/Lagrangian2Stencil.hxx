//
//  Lagrangian2Stencil.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 17/08/2018.
//

#ifndef Lagrangian2Stencil_hxx
#define Lagrangian2Stencil_hxx

template<typename T> auto StencilQuadLinLag2<T>::
GetGuess(IndexCRef index) const -> DistanceGuess {
    NormType norm = GetNorm(index);
    const ScalarType h = param.gridScale;
    norm.m *= square(h);
    norm.w *= h;
    return norm;
}

template<typename T> void StencilQuadLinLag2<T>::
SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) {
    const NormType & norm = GetNorm(index);
    l.clear();
    l.insert_after(l.before_begin(),{OffsetType(1,0),OffsetType(0,1),OffsetType(-1,0),OffsetType(0,-1),OffsetType(1,0)});
    SternBrocotRefine([&norm](OffsetCRef u, OffsetCRef v) -> bool {
        return norm.IsAcute(VectorType::CastCoordinates(u), VectorType::CastCoordinates(v));}, l);
    stencil.insert(stencil.end(),l.begin(),l.end());
    stencil.pop_back();
}

template<typename T> auto StencilQuadLinLag2<T>::
HopfLaxUpdate(IndexCRef index, const OffsetVal3 & offsetVal) -> std::pair<ScalarType,int> {
    const NormType & base = GetNorm(index);
    const ScalarType h = param.gridScale;
    
    // Try from the center offset
    assert(!offsetVal.empty());
    typedef typename NormType1::VectorType VectorType1;
    NormType1 norm1;
    const VectorType off0 = h*VectorType::CastCoordinates(offsetVal[0].first);
    norm1.m(0,0) = base.m.SquaredNorm(off0);
    norm1.w[0] = base.w.ScalarProduct(off0);
    
    const VectorType1 d1{offsetVal[0].second};
    VectorType1 grad1;
    ScalarType value = HopfLaxMinimize(norm1,d1,grad1);
    DiscreteType active = 0;
    
    if(offsetVal.size()==1) return {value,active};
    
    // Try from the first offset pair
    assert(offsetVal.size()>=2);
    NormType norm;
    norm.m(0,0) = norm1.m(0,0);
    norm.w[0] = norm1.w[0];

    VectorType off1 = h*VectorType::CastCoordinates(offsetVal[1].first);
    norm.m(0,1) = base.m.ScalarProduct(off0,off1);
    norm.m(1,1) = base.m.ScalarProduct(off1,off1);
    norm.w[1] = base.w.ScalarProduct(off1);
    
    VectorType d{d1[0],offsetVal[1].second};
    VectorType grad;
    ScalarType candidate = HopfLaxMinimize(norm,d,grad);
    
    if(candidate<value) {
        value=candidate;
        active=1;
        assert(grad.AreAllCoordinatesNonNegative());
    }
    
    if(offsetVal.size()==2) return {value,active};
    
    // try from the second offset pair
    assert(offsetVal.size()==3);
    
    off1 = h*VectorType::CastCoordinates(offsetVal[2].first);
    norm.m(0,1) = base.m.ScalarProduct(off0,off1);
    norm.m(1,1) = base.m.ScalarProduct(off1,off1);
    norm.w[1] = base.w.ScalarProduct(off1);
    d[1]=offsetVal[2].second;
    
    candidate = HopfLaxMinimize(norm,d,grad);
    if(candidate<value){
        value=candidate;
        active=2;
        assert(grad.AreAllCoordinatesNonNegative());
    }
    
    return {value,active};
}

template<typename T> auto StencilQuadLinLag2<T>::
HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow) -> RecomputeType {
    const NormType & base = GetNorm(index);
    const ScalarType h = param.gridScale;

    assert(!flow.empty());
    if(flow.size()==1){
        const VectorType off = h*VectorType::CastCoordinates(flow[0].offset);
        const ScalarType value = flow[0].weight+ base.Norm(off); // Flow[0].weight initially stores the value at neighbor
        flow[0].weight = 1;
        return {value,0.};
    }
    
    assert(flow.size()==2);
    const VectorType
    off0 = h*VectorType::CastCoordinates(flow[0].offset),
    off1 = h*VectorType::CastCoordinates(flow[1].offset);

    NormType norm;
    norm.m(0,0) = base.m.SquaredNorm(off0);
    norm.m(0,1) = base.m.ScalarProduct(off0, off1);
    norm.m(1,1) = base.m.SquaredNorm(off1);
    norm.w[0] = base.w.ScalarProduct(off0);
    norm.w[1] = base.w.ScalarProduct(off1);
    
    const VectorType d{flow[0].weight,flow[1].weight}; // Values at the neighbors
    VectorType g;
    ScalarType value = HopfLaxMinimize(norm,d,g);
    
    // If minimum is attained in the segment interior, then return.
    if(value<Traits::Infinity()){
        ScalarType width(0.);
        for(int i=0; i<Dimension; ++i){
            width += g[i]*std::fabs(value-d[i]);
            flow[i].weight = g[i];
        }
        return RecomputeType{value,width};
    }
    
    // Otherwise test the two endpoints
    const ScalarType
    val0 = norm.Norm(VectorType{1.,0.}) + d[0],
    val1 = norm.Norm(VectorType{0.,1.}) + d[1];
    if(val0<=val1){
        value=val0;
    } else {
        value=val1;
        flow[0].offset = flow[1].offset;
    }
    flow.resize(1);
    flow[0].weight = 1.;

    return RecomputeType{value,0.};
}
#endif /* Lagrangian2Stencil_hxx */
