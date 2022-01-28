//
//  EulerianStencil.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/08/2018.
//

#ifndef EulerianStencil_hxx
#define EulerianStencil_hxx

// Active stencil manipulation. iMax manually encoded in base 2

template<typename TDiff, int nSym, int nFor, int nM> int
EulerianStencil<TDiff,nSym,nFor,nM>::GetIMax(ActiveNeighFlagType b) {
    int iMax=0;
    for(int i=0; i<nMaxBits; ++i) if(b[nNeigh+i]) iMax += (1<<i);
    assert(0<=iMax && iMax<nMax);
    return iMax;
}

template<typename TDiff, int nSym, int nFor, int nM> void
EulerianStencil<TDiff,nSym,nFor,nM>::SetIMax(ActiveNeighFlagType & b, int iMax) {
    assert(0<=iMax && iMax<nMax);
    for(int i=0; i<nMaxBits; ++i) b[nNeigh+i] = iMax & (1<<i);
}

// Elementary Quadratic form


template<typename TS, size_t n> auto
QuadraticMax<TS,n>::Solve() const -> std::pair<ScalarType, int> {
    
    // Single pass, less data ?
    int iMin(0);
    ScalarType result = Infinity();
    for(int i=0; i<n; ++i){
        const DataType & d = data[i];
        if(d.a==0) continue;
        const ScalarType delta2 = d.b*d.b - d.a*d.c;
		// Mathematically, one has the guarantee that delta2>=0.
		// However this can fail due to roundoff errors, hence we consider max(0,delta2)
        const ScalarType delta = sqrt(std::max(0.,delta2));
        const ScalarType res = (d.b+delta)/d.a;
        if(n==0) return {minVal+res,0};
        if(res>=result) continue;
        iMin=i;
        result=res;
    }
    return {minVal+result,iMin};
}

template<typename TS, size_t n> void
QuadraticMax<TS,n>::
Add(ScalarType value, ScalarType weight, int i) {
    const ScalarType w=weight, v=value-minVal;
    const ScalarType wv=w*v, wv2 = (w*v)*v;
    auto & d = data[i];
    d.a+=w;
    d.b+=wv;
    d.c+=wv2;
}

template<typename TS,size_t n> void
QuadraticMax<TS,n>::Add(ScalarType value, ScalarType weight) {
    const ScalarType w=weight, v=value-minVal;
    const ScalarType wv=w*v, wv2 = (w*v)*v;
    for(auto & d : data){
        d.a+=w;
        d.b+=wv;
        d.c+=wv2;
    }
}

// Printing

template<typename TDiff, int nSym, int nFor, int nM> void
EulerianStencil<TDiff,nSym,nFor,nM>::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayRecursiveArrow(forward, 1)
    ExportArrayRecursiveArrow(symmetric, 1)
    << "}";
}

template<typename TOff, typename TScal, int VMult> std::ostream &
operator << (std::ostream & os, const EulerianDifference<TOff, TScal, VMult> & diff){
    return os << "{" << diff.baseWeight << "," << diff.offset  << "}";
}

template<typename TS, size_t n> void
QuadraticMax<TS,n>::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(minVal)
    ExportArrayArrow(data)
    << "}";
}

template <typename TS, size_t n> void
QuadraticMax<TS,n>::DataType::PrintSelf(std::ostream & os) const {
    os << "{" << a << "," << b << "," << c << "}";
}

//  ------------ Hopf-Lax update ----------

template<typename TDiff,int nSym,int nFor,int nM> auto
EulerianStencil<TDiff,nSym,nFor,nM>::
HopfLaxUpdate(OffsetType offset, ScalarType acceptedValue, const MultiplierType & mult,
              QuadType & quad, ActiveNeighFlagType & active) const -> ScalarType {
    
    if(quad.minVal==-QuadType::Infinity()){
        assert(active==0);
        quad.minVal = acceptedValue;
    }
    
    int iNeigh=0;
    for(int iMax=0; iMax<nMax; ++iMax){
        for(const auto & diff : forward[iMax]){
            if(diff.offset==offset){
                assert(!active[iNeigh]);
                active[iNeigh]=1;
                quad.Add(acceptedValue,diff.Weight(mult),iMax);
            }
            ++iNeigh;
        }
        
        for(const auto & diff : symmetric[iMax]){
            if(diff.offset==offset){
                assert(!active[iNeigh]);
                if(!active[iNeigh+1]){
                    active[iNeigh]=1;
                    quad.Add(acceptedValue,diff.Weight(mult), iMax);
                }
            }
            if(diff.offset==-offset){
                assert(!active[iNeigh+1]);
                if(!active[iNeigh]){
                    active[iNeigh+1]=1;
                    quad.Add(acceptedValue,diff.Weight(mult), iMax);
                }
            }
            iNeigh+=2;
        }
    }
    
    const auto & [value,activeQuad] = quad.Solve();
    
    /*    // Alternatively, only insert in queue if value is strictly decreased.
     ScalarType & val = values[updatedLinearIndex];
     if(result.first>=val) return;*/
    
    SetIMax(active,activeQuad);
    return value;
    
}

template<typename TDiff,int nSym,int nFor,int nM> template<typename F> auto
EulerianStencil<TDiff,nSym,nFor,nM>::
HopfLaxRecompute(const F & getVal, const MultiplierType & mult,
				 ActiveNeighFlagType active,
                 DiscreteFlowType & discreteFlow) const
-> RecomputeType {
    assert(discreteFlow.empty());
    assert(!active.none());
	std::array<int,nActiveNeigh> orderNeigh;
//    std::fill(orderNeigh.begin(),orderNeigh.end(),0);
    CappedVector<ScalarType, nActiveNeigh> weights;
    
    const int iMax = GetIMax(active);
    assert(0<=iMax && iMax<nMax);
    int iNeigh = iMax * nSingleNeigh;
    
    // We actually store the neighbor values, not weights,
	// (and the offsets) at this point.
    auto flowPush = [&getVal,&orderNeigh,&discreteFlow](const OffsetType & offset){
        int ord = 3;
        const ScalarType val = getVal(offset,ord);
        assert(ord>=1);
        orderNeigh[discreteFlow.size()] = ord;
        discreteFlow.push_back({offset,val});
    };
    
    for(const auto & diff : forward[iMax]){
        if(active[iNeigh]){
            weights.push_back(diff.Weight(mult));
            flowPush(diff.offset);
        }
        ++iNeigh;
    }
    for(const auto & diff : symmetric[iMax]){
        if(active[iNeigh]){
            assert(!active[iNeigh+1]);
            weights.push_back(diff.Weight(mult));
            flowPush(diff.offset);
        }
        if(active[iNeigh+1]){
            assert(!active[iNeigh]);
            weights.push_back(diff.Weight(mult));
            flowPush(-diff.offset);
        }
        iNeigh+=2;
    }
    
    assert(!discreteFlow.empty());
    
    QuadraticMax<ScalarType, 1> quad; quad.minVal = QuadType::Infinity();
    for(const auto & offsetValue : discreteFlow){
        quad.minVal=std::min(quad.minVal, offsetValue.weight);}
	
	// TODO : third order
    for(size_t i=0; i<discreteFlow.size(); ++i){
		const int ord = orderNeigh[i];
		assert(ord!=0);
		constexpr std::array<ScalarType,4> mult2 = {0.,1.,9./4.,121./36.};
        quad.Add(discreteFlow[i].weight,mult2[ord]*weights[i]);
    }
    
    RecomputeType result;
    result.value = quad.Solve().first;
    result.width = 0.; 
    
    for(size_t i=0; i<discreteFlow.size(); ++i){
        // Difference should already be non-negative, by construction, without sndOrder.
        const ScalarType weightPosDiff = weights[i]*std::max(0.,result.value-discreteFlow[i].weight);
        discreteFlow[i].weight = weightPosDiff;
        result.width+=weightPosDiff;
		const int ord = orderNeigh[i];
		assert(ord!=0);
		constexpr std::array<ScalarType,4> mult = {0.,1.,3./2.,11./6.};
        discreteFlow[i].weight *= mult[ord];
    }
    const ScalarType weightsSum = std::accumulate(weights.begin(), weights.end(), 0.);
    if(weightsSum>0) {result.width/=weightsSum;}
    return result;
}

#endif /* EulerianStencil_h */
