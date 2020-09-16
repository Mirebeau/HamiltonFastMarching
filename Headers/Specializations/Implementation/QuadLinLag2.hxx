//
//  QuadLinLag2.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 20/11/2018.
//

#ifndef QuadLinLag2_hxx
#define QuadLinLag2_hxx

template<typename T> template<typename Norm, bool dualize>
struct StencilQuadLinLag2<T>::MetricCaster : StencilQuadLinLag2<T>::MetricType {
	using Traits = T;
	using StencilType = StencilQuadLinLag2<T>;
	Redeclare5Types(StencilType,IndexType,MetricElementType,ScalarType,SymmetricMatrixType,NormType)
	using MetricSource = typename Traits::template DataSource<Norm>;
	std::unique_ptr<MetricSource> pMetric;
	MetricCaster(std::unique_ptr<MetricSource> pMetric_):pMetric(std::move(pMetric_)){};
	virtual bool CheckDims(const IndexType & dims) const {return pMetric->CheckDims(dims);};
	virtual MetricElementType operator()(const IndexType & index) const {
		using Sym = SymmetricMatrixType;
		using Vec = VectorType;
		if constexpr(std::is_same_v<Norm,ScalarType>) {
			ScalarType s = (*pMetric)(index);
			s = dualize ? 1./(s*s) : (s*s);
			return MetricElementType{s*Sym::Identity(),Vec::Constant(0.)};
		} else if constexpr(std::is_same_v<Norm,SymmetricMatrixType>) {
			SymmetricMatrixType s = (*pMetric)(index);
			if(dualize) {s = s.Inverse();}
			return MetricElementType{s,Vec::Constant(0.)};
		} else if constexpr(std::is_same_v<Norm,MetricElementType>) {
			const auto [m,w] = (*pMetric)(index);
			NormType norm{m,w};
			if(dualize) {norm=norm.DualNorm();}
			return MetricElementType{norm.m,norm.w};
		} else {
			static_assert(dependent_false_v<T>,"Unsupported norm type.");
		}
	};
};


template<typename T> void StencilQuadLinLag2<T>::
Setup(HFMI * that){
	Superclass::Setup(that);
	param.Setup(that);
	gridScales = SymmetricMatrixType::RankOneTensor(VectorType::FromOrigin(param.gridScales));
	
	auto & io = that->io;
	if(io.HasField("dualMetric")) {
		if(io.HasField("metric")) ExceptionMacro("Error: both primal and dual metric provided");
		SetMetricCaster<true>(that);
	} else {SetMetricCaster<false>(that);}
	cosAngleMin = io.template Get<ScalarType>("cosAngleMin",cosAngleMin);
	if(io.template Get<ScalarType>("refineStencilAtWallBoundary",0.) && io.HasField("walls")){
		wallBoundaryAngularResolution = io.template Get<ScalarType>("wallBoundaryAngularResolution",wallBoundaryAngularResolution,2);
		const DomainType & dom = that->pFM->dom;
		pDom = &dom;
		walls.dims = that->pFM->stencilData.dims;
		walls.resize(walls.dims.Product());
		auto pWalls = that->template GetIntegralField<bool>("walls");
		for(DiscreteType i=0; i<walls.size(); ++i){walls[i]=(*pWalls)(walls.Convert(i));}
	}
}

template<typename T> template<bool dualize>
void StencilQuadLinLag2<T>::SetMetricCaster(HFMI * that){
	using Sym = SymmetricMatrixType;
	using Met = MetricElementType;
	const std::string key = dualize ? "dualMetric" : "metric";
	const int elemSize = that->FieldElementSize(key);
	switch(elemSize){
		case 1:
			pMetric.reset(new MetricCaster<ScalarType,dualize>(that->template GetField<ScalarType>(key,false)));
			break;
		case Sym::InternalDimension:
			pMetric.reset(new MetricCaster<Sym,dualize>(that->template GetField<Sym>(key,false)));
			break;
		default:
			if(dualize) {pMetric.reset(new MetricCaster<Met,dualize>(that->template GetField<Met>(key,false)));}
			else {pMetric = that->template GetField<Met>(key,false);}
	}

}

template<typename T> auto StencilQuadLinLag2<T>::
GetNorm(IndexCRef index) const -> NormType {
	assert(pMetric!=nullptr);
	return Rescale((*pMetric)(index));
}

template<typename T> auto StencilQuadLinLag2<T>::
GetGuess(const PointType & p) const -> DistanceGuess {
	assert(pMetric!=nullptr);
	return Rescale(MapWeightedSum<MetricElementType>(*pMetric,this->pFM->dom.Neighbors(p)));
}

template<typename T> auto StencilQuadLinLag2<T>::
Rescale(const MetricElementType & data)
const -> NormType {
	const auto & [m,v] = data;
	NormType norm{m,v};
	
	norm.m = norm.m.ComponentWiseProduct(gridScales);
	norm.w = norm.w.ComponentWiseProduct(VectorType::FromOrigin(param.gridScales));
	return norm;
}


template<typename T> void StencilQuadLinLag2<T>::
SetNeighbors(IndexCRef index, std::vector<OffsetType> & stencil) {
	const NormType & norm = GetNorm(index);
	assert(tmp_stencil.empty());
	tmp_stencil.insert(tmp_stencil.end(),{OffsetType(1,0),OffsetType(0,-1),OffsetType(-1,0),OffsetType(0,1)});

	if(OnWallBoundary(index)){
		auto pred = [&norm,this,&index](OffsetCRef u, OffsetCRef v) -> bool {
			const VectorType vu = VectorType::CastCoordinates(u), vv = VectorType::CastCoordinates(v);
			IndexType pu = index + IndexDiff::CastCoordinates(u), pv = index + IndexDiff::CastCoordinates(v);
			if(!this->pDom->Periodize(pu,index).IsValid() || !this->pDom->Periodize(pv,index).IsValid()) return true;
			const bool wu = this->walls(pu), wv = this->walls(pv);
			if(wu && wv) return true;
			if(!wu && !wv) return CosAngle(norm, vu, vv)>=this->cosAngleMin;
			return square(LinearAlgebra::Determinant(vu,vv)/this->wallBoundaryAngularResolution) > vu.SquaredNorm()*vv.SquaredNorm();
		};
		
		SternBrocotRefine(pred, stencil, tmp_stencil);
	} else {
		/* // Predicate based version. Uses two or three gradient evaluations per point.
		auto pred = [&norm,this](OffsetCRef u, OffsetCRef v) -> bool {
			return CosAngle(norm, VectorType::CastCoordinates(u),
							VectorType::CastCoordinates(v))>=this->cosAngleMin;};
		
		SternBrocotRefine(pred, stencil, tmp_stencil);
		 */
		
		SternBrocotRefine_AcuteBound(norm, cosAngleMin, stencil,
									 tmp_stencil, tmp_stencil_vec, tmp_stencil_scal);

	}
	
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
	
	// Try from the center offset
	assert(!offsetVal.empty());
	typedef typename NormType1::VectorType VectorType1;
	NormType1 norm1;
	const VectorType off0 = VectorType::CastCoordinates(offsetVal[0].first);
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
	
	VectorType off1 = VectorType::CastCoordinates(offsetVal[1].first);
	norm.m(0,1) = base.m.ScalarProduct(off0,off1);
	norm.m(1,1) = base.m.ScalarProduct(off1,off1);
	norm.w[1] = base.w.ScalarProduct(off1);
	
	VectorType d{d1[0],offsetVal[1].second};
	VectorType grad;
	ScalarType candidate = HopfLaxMinimize(norm,d,grad);
	
	if(candidate<value) {
		value=candidate;
		active=1;
		assert(grad.IsNonNegative());
	}
	
	if(offsetVal.size()==2) return {value,active};
	
	// try from the second offset pair
	assert(offsetVal.size()==3);
	
	off1 = VectorType::CastCoordinates(offsetVal[2].first);
	norm.m(0,1) = base.m.ScalarProduct(off0,off1);
	norm.m(1,1) = base.m.ScalarProduct(off1,off1);
	norm.w[1] = base.w.ScalarProduct(off1);
	d[1]=offsetVal[2].second;
	
	candidate = HopfLaxMinimize(norm,d,grad);
	if(candidate<value){
		value=candidate;
		active=2;
		assert(grad.IsNonNegative());
	}
	
	return {value,active};
}

template<typename T> auto StencilQuadLinLag2<T>::
HopfLaxRecompute(IndexCRef index, DiscreteFlowType & flow) -> RecomputeType {
	const NormType & base = GetNorm(index);
	
	assert(!flow.empty());
	if(flow.size()==1){
		const VectorType off = VectorType::CastCoordinates(flow[0].offset);
		const ScalarType offsetNorm = base.Norm(off);
		const ScalarType value = flow[0].weight+ offsetNorm; // Flow[0].weight initially stores the value at neighbor
		flow[0].weight = 1/offsetNorm;
		return {value,0.};
	}
	
	assert(flow.size()==2);
	const VectorType
	off0 = VectorType::CastCoordinates(flow[0].offset),
	off1 = VectorType::CastCoordinates(flow[1].offset);
	
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
		assert(g.Sum()>0);
		return RecomputeType{value,width/g.Sum()};
	}
	
	// Otherwise test the two endpoints
	const ScalarType
	norm0 = norm.Norm(VectorType{1.,0.}),
	norm1 = norm.Norm(VectorType{0.,1.});
	const ScalarType
	val0 = norm0 + d[0],
	val1 = norm1 + d[1];
	if(val0<=val1){
		value=val0;
		flow[0].weight = 1./norm0;
	} else {
		value=val1;
		flow[0].offset = flow[1].offset;
		flow[0].weight = 1./norm1;
	}
	flow.resize(1);

	return RecomputeType{value,0.};
}

#endif /* QuadLinLag2_h */
