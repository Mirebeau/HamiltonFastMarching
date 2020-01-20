#ifndef SeismicNorm_hpp
#define SeismicNorm_hpp

// *********** Array access ************
template<typename TC, size_t VD> int 
SeismicNorm<TC, VD>::VoigtIndex(int i, int j) const {
	[[maybe_unused]] size_t const d = Dimension;
	assert(0 <= i && i < d && 0 <= j && j < d);
	switch (Dimension) {
	case 1: return 0;
	case 2: {
		static const std::array<std::array<int, 2>, 2>
		indices = { { {{0,2}},{{2,1}} } };
		return indices[i][j];
	}
	case 3: {
		static const std::array<std::array<int, 3>, 3>
		indices = { { {{0,5,4}},{{5,1,3}}, {{4,3,2}} } };
		return indices[i][j];
		}
	default: assert(false); // dimension non supportée
		return std::numeric_limits<int>::max();
	}
};

template<typename TC, size_t VD> std::pair<int,int>
SeismicNorm<TC, VD>::VoigtIndexInv(int k) const {
	using int2 = std::pair<int,int>;
	[[maybe_unused]] size_t const d = Dimension;
	assert(0 <= k && k < d*(d + 1) / 2);
	switch (Dimension) {
		case 1: return { 0,0 };
		case 2: {
			static const std::array<int2, 3>
			indices = {{ {0,0}, {1,1}, {0,1} }};
			return indices[k];
		}
		case 3: {
			static const std::array<int2,6>
			indices = {{ {0,0}, {1,1}, {2,2}, {1,2}, {0,2}, {0,1} }};
			return indices[k];
		}
		default: {
			assert(false); // unsupported Dimension
			static const int maxnum=std::numeric_limits<int>::max();
			return { maxnum, maxnum };
		}
	}// switch Dimension
};

//Constructor
template<typename TC, size_t VD>
SeismicNorm<TC,VD>::SeismicNorm(SymmetricMatrixType const& m) {
	for (int a = 0; a < VoigtDimension; ++a) {
		for (int b = a; b < VoigtDimension; ++b) {
			const auto [i, j] = VoigtIndexInv(a);
			const auto [k, l] = VoigtIndexInv(b);
			hookeTensor(a, b) = m(i, j)*m(k, l);
		}
	}
}

// Tensor contraction

template<typename TC, size_t VD>
template<typename T, typename TOut> auto
SeismicNorm<TC, VD>::Contract(Vector<T, Dimension> const& u)
const -> SymmetricMatrix<TOut, Dimension> {
	int const d = Dimension;
	auto m = SymmetricMatrix<TOut, d>::Zero();
	for (int i = 0; i < d; ++i) {
		for (int j = 0; j < d; ++j) {
			for (int k = i; k < d; ++k) {
				for (int l = 0; l < d; ++l) {
					m(i, k) += Coefficient(i, j, k, l)*u[j] * u[l];
				}
			}
		}
	}
	return m;
}


template<typename TC, size_t VD>
template<typename T, typename TOut, typename S> TOut
SeismicNorm<TC, VD>::Constraint(Vector<T, Dimension> const& u, S const & s)
const {
	int const d = Dimension;
	typedef decltype(u[0]*u[0]) T2;
	const auto m = Contract<T,T2>(u);
	typedef decltype(u[0]*u[0]-s) T2b;
	Vector<T2b,Dimension> diag;
	for(int i=0; i<d; ++i) diag[i]=m(i,i)-s;
	static_assert(Dimension <= 3, "Error : unsupported dimension");
	
	// Determinant of opposite
	if constexpr (Dimension == 1) {
		return - diag[0];
	} else if constexpr (Dimension == 2) {
		return diag[0]*diag[1] - m(1, 0)*m(1, 0); //m(0, 0)*m(1, 1)
	}
	else if constexpr (Dimension == 3) {
//		return
//		-(diag[0]*diag[1]*diag[2] + m(0,1)*m(1,2)*m(2,0) + m(0,2)*m(1,0)*m(2,1)
//		- m(0,2)*diag[1]*m(2,0) - m(0,1)*m(1,0)*diag[2] - diag[0]*m(1,2)*m(2,1));
		
		return  // Optimized expression
		diag[0]*(m(1,2)*m(2,1)-diag[1]*diag[2])
		+m(0,1)*(m(1,0)*diag[2]-ComponentType(2.)*m(1,2)*m(2,0))
		+m(0,2)*(diag[1]*m(2,0));
	};
//	return (d%2==0 ? m.Determinant() : -m.Determinant());
}

// ********* Acceleration of GradNorm computation *******

template<typename TC, size_t VD> auto
SeismicNorm<TC, VD>::AD2Constraint_tmp() const
-> AD2SymType {
	// Optimized version
	AD2SymType s;
	for(int i=0; i<Dimension; ++i){
		for(int k=0; k<=i; ++k){
			for(int j=0; j<Dimension; ++j){
				for(int l=0; l<=j; ++l){
					s(i,k).m(j,l) = Coefficient(i, j, k, l) + Coefficient(i,l,k,j);
					// Optimization opportunity : in dimension d=2, the above two are identical unless
					//(i==1 && j==1 && k==0 && l==0)
				}
			}
		}
	}
	return s;
}

template<typename TC, size_t VD> auto
SeismicNorm<TC, VD>::AD2Constraint_zero(AD2SymType const & s) const
-> AD2Type {
	AD2Type c;
	c.x=ComponentType(1.);
	c.v.fill(ComponentType(0.));
	
	// c.m = -sum_i s(i,i).m
	c.m = s(0,0).m;
	for(int i=0; i<Dimension; ++i) {c.m += s(i,i).m;}
	c.m = -c.m;
	return c;
}

template<typename TC, size_t VD> auto
SeismicNorm<TC, VD>::AD2Constraint(VectorType const & p, AD2SymType & s) const
-> AD2Type {
	AD2Type c;
	for(int i=0; i<Dimension; ++i){
		for(int j=0; j<=i; ++j){
			s(i,j).v = s(i,j).m*p;
			s(i,j).x = 0.5*(s(i,j).v.ScalarProduct(p));
		}
	}
	for(int i=0; i<Dimension; ++i) {s(i,i).x -= 1.;}
	
	// Determinant of opposite
	if constexpr(Dimension==1) {
		c = -s(0,0);
	} else if constexpr(Dimension==2){
		c = s(0,0)*s(1,1)-s(1,0)*s(1,0);
	} else if constexpr(Dimension==3) {
		
		//c= (s(0,2)*s(1,1)*s(2,0) + s(0,1)*s(1,0)*s(2,2) + s(0,0)*s(1,2)*s(2,1))
		//-(s(0,0)*s(1,1)*s(2,2) + s(0,1)*s(1,2)*s(2,0) + s(0,2)*s(1,0)*s(2,1));
	
		 //c=
		 //s(0,0)*(s(1,2)*s(2,1)-s(1,1)*s(2,2))
		 //+s(0,1)*(s(1,0)*s(2,2)-s(1,2)*s(2,0))
		 //+s(0,2)*(s(1,1)*s(2,0)-s(1,0)*s(2,1));
		
		c= // Optimized expression
		s(0,0)*(s(1,2)*s(2,1)-s(1,1)*s(2,2))
		+s(0,1)*(s(1,0)*s(2,2)-ComponentType(2.)*s(1,2)*s(2,0))
		+s(0,2)*(s(1,1)*s(2,0));

	} else {
		assert(false);
	}
	return c;
	
	/*
	 // For reference, the following code based on (second order, forward) automatic differentiation.
	 // Produced mathematically identical results. However, it is a bit slow.
	 
	 
	 using AD2VectorType = Vector<AD2Type, Dimension>;
	 AD2VectorType pAD;
	 for (int j = 0; j < Dimension; ++j) {
	 pAD[j] = AD2Type(p[j], j);
	 }
	 return Constraint(pAD);
	 
	 */
}

// ****** Primal norm computation *******
// The vast majority of CPU time goes here.
// This procedure has been shown to fail if anellipticity is too pronounced.
// In that case, a path following method would be appropriate.
// Seems to work in realistic cases, though. (Failure observed with anell = 0.0011..., 2D)
template<typename TC, size_t VD> auto
SeismicNorm<TC, VD>::GradNorm(VectorType const& q) const
-> VectorType {
	constexpr bool isSimple = std::is_same_v<decltype(q[0]==q[0]),bool>;
	if constexpr(isSimple) {assert(q.SquaredNorm()!=0);}
	using std::sqrt;
	VectorType p = VectorType::Constant(ComponentType(0.));
	
	// Perform one step of sequential quadratic programming
	// Aim for maximizing <q,p> subject to constraint >=0, differentiated at p
	auto sqp_step = [&p,&q](const AD2Type & c){
		const SymmetricMatrixType d = c.m.Inverse();
		const VectorType dv = d*c.v, dq = d*q;
		const ComponentType num = dv.ScalarProduct(c.v) - 2.*c.x;
		const ComponentType	den = dq.ScalarProduct(q);
//		const ComponentType num = d.SquaredNorm(c.v) - 2.*c.x;
//		const ComponentType	den = d.SquaredNorm(q);
		
		// MSVC requires to redeclare isSimple
		constexpr bool isSimple = std::is_same_v<decltype(q[0] == q[0]), bool>;
		if constexpr(isSimple) {assert(num*den >= 0);}
		const ComponentType lambda = -sqrt(num / den);
		// Optimization (??) : saved d*c.v and d*q from the squared norm computation,
		// and avoid one matrix-vector product below
		const VectorType h = lambda*dq - dv;
		//const VectorType h = d * (lambda*q - c.v);
		p += h;
	};
	
	AD2SymType constraint_tmp = AD2Constraint_tmp();
	for(int n=0; n<nIterMax_GradNorm; ++n){
		if constexpr(isSimple) {assert(n>0 || p.IsNull());}
		AD2Type c =
		n==0 ?  AD2Constraint_zero(constraint_tmp) :
		(useMetilOptimization && Dimension==3) ? AD2Constraint_Metil(p) :
		AD2Constraint(p, constraint_tmp);
		
		sqp_step(c);
		
		// Optimization : adaptive stopping criterion
		using std::abs;
		if constexpr(isSimple) {if(abs(c.x) < constraintBound_GradNorm) break;}
		else {
#ifdef XSIMD_HPP
			if( xsimd::all(abs(c.x) < constraintBound_GradNorm) ) break;
#else
			static_assert(isSimple);
#endif
		}
	}

	return p;
}

// ************* Dual norm computation **********
template<typename TC, size_t VD>
template<typename T>
T SeismicNorm<TC, VD>::DualNorm(Vector<T,Dimension> const & u) const {
	using Pol1_1 = Pol1<T,1>;
	using Pol1_d = Pol1<T,Dimension>;
	Pol1_1 lambda;
	lambda.c[0] = T(0.);
	lambda.c[1] = T(1.);

	// Get Characteristic polynomial
	Pol1_d const chi = Constraint<T,Pol1_d,Pol1_1>(u,lambda);
	
	// Solve
	using std::sqrt;
	T norm2; // Largest root of chi
	const auto & c = chi.c;
	
	if constexpr(Dimension==1){
		norm2 = -c[0]/c[1];
	}
	
	if constexpr(Dimension==2) {
		const T b = c[1]/2.;
		const T delta = b*b - c[0]*c[2];
		assert(delta>=0.);
		norm2 = (-b+sqrt(delta))/c[2];
	}
	
	if constexpr(Dimension==3){
		assert(false);
	};
	assert(norm2>=0.);
	return sqrt(norm2);
}

template<typename TC, size_t VD> auto
SeismicNorm<TC, VD>::GradDualNorm(VectorType const & u)
const -> VectorType {
	using AD1Type = DifferentiationType<ComponentType, VectorType>;
	Vector<AD1Type, Dimension> adU;
	for(int i=0; i<Dimension; ++i){
		adU[i] = AD1Type(u[i],i);}
	
	return DualNorm(adU).v;
}


// ------ Anellipticity -----------

template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
EllipticApproximation() const -> SymmetricMatrixType {
	auto res = SymmetricMatrixType::Zero();
	for(int i=0; i<Dimension; ++i){
		for(int j=i; j<Dimension; ++j){
			for(int k=0; k<Dimension; ++k){
				res(i,j) += Coefficient(i, k, j, k);
			}
		}
	}
	return res;
}


// Essentially, a measure of non-smoothness
template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
Anellipticity() const -> ComponentType {
	static_assert(Dimension==2);
	
	using std::sqrt, std::min;
	using AffType = Pol1<ComponentType, 1>;
	using QuadType = Pol1<ComponentType, 2>;
	Vector<AffType, Dimension> u;
	u[0].c = {1.,0.};
	u[1].c = {0.,1.};
	
	// Matrix is full of second order polynomials
	const auto s = Contract<AffType,QuadType>(u);
	const auto
	num1 = (s(0,0)-s(1,1))/2.,
	& num2 = s(0,1),
	den = (s(0,0)+s(1,1))/2.;
	
	auto anell = [&num1, &num2, &den](ComponentType const r){
		ComponentType const n1=num1(r),n2=num2(r),d=den(r);
		assert(d>0);
		return sqrt(n1*n1+n2*n2)/d;
	};
	// A good anellipticity measure is the inf of sqrt(num1^2+num2^2)/den.
	// Problems arise when this quantity vanishes.
	// Here we content ourselves with a basic upper bound, solving for num1=0
	// TODO : is there any guarantee about this ?

	auto anellEval = [&anell](QuadType const & num){
		ComponentType const c = num.c[0], b=num.c[1]/2, a=num.c[2];
		assert(a!=0);
		// TODO: handle this case.
		ComponentType const delta = b*b-a*c;
		if(delta>=0){
			ComponentType const sdelta = sqrt(delta);
			ComponentType const r0 = (-b-sdelta)/a, r1=(-b+sdelta)/a;
			return min(anell(r0),anell(r1));
		} else {
			ComponentType const r=-b/a;
			return anell(r);
		}
	};

	return min(anellEval(num1),anellEval(num2));
}

#endif // SeismicNorm_hpp
