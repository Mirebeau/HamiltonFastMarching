//
//  TTINorm.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 02/01/2020.
//

#ifndef TTINorm_hxx
#define TTINorm_hxx

// ----- Printing -----

template<typename TS,int VD> void TTINorm<TS,VD>::
PrintSelf(std::ostream & os) const {
	os << "{"
	ExportVarArrow(linear)
	ExportVarArrow(quadratic)
	ExportVarArrow(transform)
	<< "}";
}

template<typename TS,int VD> void TTINorm<TS,VD>::Properties::
PrintSelf(std::ostream & os) const {
	os << "{"
	ExportVarArrow(optimDirection)
	ExportVarArrow(tMin)
	ExportVarArrow(tMax)
	<< "}";
}

// -------------

template<typename TS, int VD> template<typename T>
T TTINorm<TS,VD>::Level(const Vector<T,Dimension> & p0) const {
	using Vec = Vector<T,Dimension>;
	Vec p = transform.Inverse().Transpose()*p0;
	for(auto & x : p) x*=x;
	if constexpr(Dimension==2){
		return 1.-(p.ScalarProduct(linear) + 0.5*quadratic.SquaredNorm(p));
	} else { static_assert(Dimension==3,"Unsupported dimension");
		Vector<T,2> q{p[0]+p[1],p[2]};
		return 1.-(q.ScalarProduct(linear) + 0.5*quadratic.SquaredNorm(q));
	}
}

template<typename TS, int VD> auto
TTINorm<TS,VD>::smallest_positive_root(ScalarType l, ScalarType q) const -> ScalarType {
	// Returns smallest positive s such that 0.5*q*s^2 + l*s -1 = 0
	// Note that for physical TTI coefficients, both roots are positive, see the paper.
	if(q==0){assert(l>0); return 1/l;}
	const ScalarType delta = l*l + 2.*q;
	assert(delta>=0);
	using std::sqrt;
	const ScalarType sdelta = sqrt(delta);
	const ScalarType rm = (- l - sdelta)/q, rp = (- l + sdelta)/q;
	const ScalarType rmin=std::min(rm,rp), rmax=std::max(rm,rp);
	assert(rmax>0);
	return rmin>0 ? rmin : rmax;
}
									   
template<typename TS, int VD> auto
TTINorm<TS,VD>::Props() const -> Properties {
	using std::sqrt;
	Properties result;
	auto slope = [](const LinearType & l) -> ScalarType {
		// returns t such that l is proportionnal to (1-t,t)
		assert(l.IsNonNegative());
		return l[1]/l.Sum();
	};
	
	// Solve along axis 0
	// equation
	const ScalarType root0 = smallest_positive_root(linear[0],quadratic(0,0));
	const ScalarType root1 = smallest_positive_root(linear[1],quadratic(1,1));
	
	result.tMin = slope({linear[0]+ quadratic(0,0)*root0,linear[1]+quadratic(1,0)*root0});
	result.tMax = slope({linear[0]+ quadratic(0,1)*root1,linear[1]+quadratic(1,1)*root1});
	
	if(result.tMin < result.tMax){
		result.optimDirection = -1;
	} else {
		std::swap(result.tMin,result.tMax);
		result.optimDirection=1;
	}
	
	return result;
}

template<typename TS, int VD> auto TTINorm<TS,VD>::
Selling(ScalarType t) const -> SellingPath {
	// Generate the extremal symmetric matrices
	const TransformType & A = transform.Inverse();
	using Sym = SymmetricMatrixType;
	const Sym D0 = Dimension==2 ? Sym::RankOneTensor(A.Row(0)) :
	Sym::RankOneTensor(A.Row(0)) + Sym::RankOneTensor(A.Row(1));
	const Sym D1 = Sym::RankOneTensor(A.Row(Dimension-1));
	return SellingPath(D0,D1,t);
}

template<typename TS,int VD> template<typename T> T TTINorm<TS,VD>::
UpdateValue(const T & t, const NeighborValuesType & val,
			const SellingPath & path) const {
	// Sort the values
	std::array<ScalarType,SymDimension> indices;
	for(int i=0; i<SymDimension; ++i) {indices[i]=i;}
	std::sort(indices.begin(),indices.end(),
			  [&val](int i,int j){return val[i]<val[j];});
	const ScalarType valMin=val[indices[0]];
	const ScalarType inf = std::numeric_limits<ScalarType>::infinity();
	if(valMin==inf) return T(inf);
	
	// Get the multiplier
	// TODO : optimization opportunity in Multiplier(t)
	const T & mult = Multiplier(t);
	
	// Compute the update, solving ax^2-2bx+c=0
	using std::sqrt; using std::max;
	ScalarType a(0),b(0),c(-RemoveAD(mult)),sol(inf);
	int r=0;
	for(; r<SymDimension; ++r){
		const int i = indices[r];
		// Optimization opportinity : first loop yiels v=0 (but w!=0)
		const ScalarType v = val[i] - valMin;
		if(v>=sol) break;
		
		const ScalarType
		w_ = path.weights0[i]+RemoveAD(t)*(path.weights1[i]-path.weights0[i]),
		w = max(ScalarType(0), w_),
		wv=w*v, wvv=wv*v;
		
		// Mathematically, this weight is guaranteed to be non-negative.
		// However, this can fail due to roundoff errors, hence the max(0,...)
		assert(w_>=-1e-10);
		
		a+=w;
		b+=wv;
		c+=wvv;
		
		if(a <= RemoveAD(mult)*1e-10){continue;}
		const ScalarType delta = b*b-a*c;
		// Mathematically, one has the guarantee that delta>=0.
		// However, this can fail due to roundoff errors, hence the max(0,...)
		assert(delta>=-1e-10);
		const ScalarType sdelta = sqrt(max(ScalarType(0),delta));
		sol = (b+sdelta)/a;
	}
	
	// Recompute the update, if the required type is not the scalar type
	if constexpr(std::is_same_v<ScalarType,T>){return valMin+sol;}
	else {
		if(sol==inf) {return T(inf);}
		T a(0),b(0),c(-mult);
		for(int r_=0; r_<r; ++r_){
			const int i = indices[r_];
			const ScalarType v = val[i] - valMin;
			const T w = path.weights0[i]+t*(path.weights1[i]-path.weights0[i]),
			wv=w*v, wvv=wv*v;
			a+=w;
			b+=wv;
			c+=wvv;
		}
		const T delta = b*b-a*c;
		assert(delta>=0.);
		const T sol = (b+sqrt(delta))/a;
		return valMin+sol;
	} /*else {
		// optimized variant, assuming t is a canonical one dimensional AD type.
		if(sol==inf) {return T(inf);}
		ScalarType ap(0),bp(0),cp(0);
		for(int r_=1; r-<r; ++r_){
			const int i = indices[r_];
			const ScalarType v = val[i] - valMin;
			const T wp = path.weights1[i]-path.weights0[i],
			wv=w*v, wvv=wv*v;
			ap+=w;
			bp+=wv;
			cp+=wvv;
		}
		T a_,b_,c_,delta;
		using AD2Type = AD2<ScalarType,1>;
		using AD1Type = DifferentiationType<ScalarType,Vector<ScalarType, 1> >;
		if constexpr(std::is_same<T,AD2Type>){
			a_ = AD2Type{a,{ap}};
			b_ = AD2Type(b,{bp}};
						 c_ = AD2Type(c,{
		
		} else if constexpr(std::is_same<T,AD1Type>{
			a_=
		}

		
		
	}*/
}

template<typename TS, int VD> template<typename T> T TTINorm<TS,VD>::
Multiplier(const T & t) const {
	// Returns beta such that diag(1-t,t)/beta is a tangent ellipse.
	// Document : EnvelopeSeismicEqn.tex
	using Sym2 = QuadraticType;
	using Vec2 = LinearType;
	using DVec2 = Vector<T,2>;
	
	if(quadratic.data.IsNull()){ // Riemannian case
		return T(1./linear.Sum());}
	
	const Sym2 Q = quadratic.Comatrix();
	const Vec2 & l = linear;
	const Vec2 Ql = Q*l;
	const ScalarType detQ = Q.Determinant();
	const ScalarType lQl = l.ScalarProduct(Ql);
	
	constexpr bool optimized = true;
	using AD2Type = AD2<ScalarType,1>;
	if constexpr(optimized && std::is_same_v<T,AD2Type>){
		// Bypass AD2 implementation for efficiency
		assert(t.v[0]==1 && t.m(0,0)==0);
		const ScalarType t_=t.x,t=t_,s=1.-t;
		const ScalarType &
		a=Q(0,0),b=Q(0,1),c=Q(1,1),
		l0=l[0],l1=l[1];
		const ScalarType
		as=a*s, ct=c*t,
		l01=l0+l1, l1s=l1*s, l0t=l0*t, l0tl1s=l0t-l1s;
		
		const T
		vQv{as*s+2*b*s*t+ct*t,(-2.)*(as+b*(t-s)-ct),2.*(a-2.*b+c)},
		detVL2{square(l0tl1s),2.*l01*l0tl1s,2.*square(l01)},
		lQv{Ql[0]*s+Ql[1]*t,Ql[1]-Ql[0],0.};
		
		const T num=detVL2+2*vQv;
		const int signNum = num>0 ? 1 : -1;
		const T sdelta = sqrt(vQv*(2*detQ+lQl));
		const T den = signNum * sdelta + lQv;
		return num/den;
	} else {
		const DVec2 v{1-t,t};
		const DVec2 Qv = Q*v;
		const T detVL = Determinant(v, l);
		
		const T lQv = l.ScalarProduct(Qv);
		const T vQv = v.ScalarProduct(Qv);
		
		const T num = detVL*detVL + 2*vQv;
		const int signNum = num>0 ? 1 : -1;
		
		using std::sqrt;
		const T delta = vQv*(2*detQ+lQl);

		const ScalarType delta0 = RemoveAD(delta);
		if(delta0<1e-10){
			// delta is mathematically guaranteed to be non-negative, up to roundoff errors
			assert(delta0>-1e-10);
			// delta=0 is a degenerate case, where the conic is a union of two lines,
			// either parallel or intersecting (outside ot [0,inf[^2 in the latter case).
			// In that case, both factors vQv and 2*detQ+lQl vanish.
			const ScalarType root0 = smallest_positive_root(linear[0],quadratic(0,0));
			const ScalarType beta = root0*RemoveAD(v[0]);
			assert(std::abs(beta - // other way to compute same multiplier
				smallest_positive_root(linear[1],quadratic(1,1))*RemoveAD(v[1]))<1e-7);
			return T(beta);
		}
		const T sdelta = sqrt(delta);
		const T den = signNum * sdelta + lQv;
		return num/den;
	}
}

// ----------------- Norm computation --------------
template<typename TS,int VD> auto TTINorm<TS,VD>::
Gradient(const VectorType & q0) const -> VectorType {
	// Use sequential quadratic programming to solve the optimization problem
	// sup <v,w> subject to Level(w)>=0,
	// where one is implicitly restricted to the connected component of the origin.
	
	// TODO : some optimizations are possible, if this becomes a limiting factor,
	// especially in 3D (exploit transversal isotropy, to reduce to 2D).
	constexpr bool optimized = true;
	constexpr int nIter_SQP = 8;
		
	if(!optimized){
		using Diff2 = AD2<ScalarType,Dimension>;
		VectorType p = VectorType::Constant(0); // Initial guess
		for(int iter=0; iter<nIter_SQP; ++iter){
			p+=Level(Diff2::Perturbation(p)).SQP(q0);}
		return p;
	}
	
	// Optimized variant
	using Diff2 = AD2<ScalarType,2>; // Two dimensional perturbation, vector, matrix
	using Vec2 = Vector<ScalarType,2>;
	using Sym2 = SymmetricMatrix<ScalarType,2>;
	
	const ScalarType &
	a=linear[0],b=linear[1],
	c=quadratic(0,0),d=quadratic(0,1),e=quadratic(1,1);
	
	const VectorType q1 = transform*q0;
	using std::sqrt;
	const Vec2 q{Dimension==2 ? q1[0] : sqrt(q1[0]*q1[0]+q1[1]*q1[1]), q1[Dimension-1]};
	// Solve analytically the first sqp step
	const ScalarType q0a=q[0]/a,q1b=q[1]/b;
	using std::sqrt;
	Vec2 p = Vec2{q0a,q1b}/sqrt(q[0]*q0a+q[1]*q1b);
		
	for(int i=1; i<nIter_SQP; ++i){
		// Evaluate constraint and derivatives
		const ScalarType & x=p[0],y=p[1];
		const ScalarType x2=x*x,y2=y*y,
		cx2=c*x2,ey2=e*y2,dx2=d*x2,dy2=d*y2;
		const Diff2 lvl
		= Diff2{1.-(a*x2+b*y2+0.5*(cx2*x2+2.*dx2*y2+ey2*y2)),
			(-2.)*Vec2{x*(a+cx2+dy2),y*(b+dx2+ey2)},
			(-2.)*Sym2{a+3*cx2+dy2,2.*d*x*y,b+dx2+3.*ey2}
		};
		p+=lvl.SQP(q);
		
	}
	const TransformType ttrans = transform.Transpose();
	
	if constexpr(Dimension==2){return ttrans*p;}
	else {
		const ScalarType n01 = q[0]; // Norm of components 1 and 2 of q1
		if(n01==0.){return ttrans*VectorType{0.,0.,p[1]};}
		else {return ttrans*VectorType{p[0]*q1[0]/n01,p[0]*q1[1]/n01,p[1]};}
	}
	
	
}

template<typename TS,int VD> auto TTINorm<TS,VD>::
Norm(const VectorType & v) const -> ScalarType {
	// Use Euler's identity to compute the norm
	if(v.IsNull()) {return 0.;}
	else {return Gradient(v).ScalarProduct(v);}
}

#endif /* TTINorm_h */
