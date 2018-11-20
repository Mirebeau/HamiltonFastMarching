//
//  QuadLinLag2.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef QuadLinLag2_hxx
#define QuadLinLag2_hxx


template<typename T> void StencilQuadLinLag2<T>::
Setup(HFMI * that){
	Superclass::Setup(that); param.Setup(that);
	pMetric = that->template GetField<MetricElementType>("metric",false);
	auto & io = that->io;
	if(io.template Get<ScalarType>("refineStencilAtWallBoundary",0.) && io.HasField("walls")){
		wallBoundaryAngularResolution = io.template Get<ScalarType>("wallBoundaryAngularResolution",wallBoundaryAngularResolution,2);
		const DomainType & dom = that->pFM->dom;
		pDom = &dom;
		walls.dims = that->pFM->stencilData.dims;
		walls.resize(walls.dims.ProductOfCoordinates());
		auto pWalls = that->template GetIntegralField<bool>("walls");
		for(DiscreteType i=0; i<walls.size(); ++i){walls[i]=(*pWalls)(walls.Convert(i));}
	}
}

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
	
	if(OnWallBoundary(index)){
		SternBrocotRefine([&norm,this,&index](OffsetCRef u, OffsetCRef v) -> bool {
			const VectorType vu = VectorType::CastCoordinates(u), vv = VectorType::CastCoordinates(v);
			IndexType pu = index + IndexDiff::CastCoordinates(u), pv = index + IndexDiff::CastCoordinates(v);
			if(!this->pDom->Periodize(pu,index).IsValid() || !this->pDom->Periodize(pv,index).IsValid()) return true;
			const bool wu = this->walls(pu), wv = this->walls(pv);
			if(wu && wv) return true;
			if(!wu && !wv) return norm.IsAcute(vu, vv);
			return square(LinearAlgebra::Determinant(vu,vv)/this->wallBoundaryAngularResolution) > vu.SquaredNorm()*vv.SquaredNorm();
		},l);
	} else {
		SternBrocotRefine([&norm](OffsetCRef u, OffsetCRef v) -> bool {
			return norm.IsAcute(VectorType::CastCoordinates(u), VectorType::CastCoordinates(v));}, l);
	}
	stencil.insert(stencil.end(),l.begin(),l.end());
	stencil.pop_back();
	
	// TODO : add criterion for walls, etc.
	// Add passthrough in WallObstruction test
}

template<typename T> bool StencilQuadLinLag2<T>::
OnWallBoundary(IndexCRef index) const {
	if(walls.empty() || walls(index)) return false;
	for(int i=0; i<Dimension; ++i){
		for(int eps=-1; eps<=1; eps+=2){
			IndexType neigh=index;
			if(!pDom->Periodize(neigh,index).IsValid()) continue;
			if(walls(neigh)) return true;
		}
	}
	return false;
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

#endif /* QuadLinLag2_h */
