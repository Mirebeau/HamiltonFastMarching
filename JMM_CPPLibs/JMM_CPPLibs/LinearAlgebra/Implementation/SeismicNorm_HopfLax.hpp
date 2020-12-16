//
//  SeismicNorm_HopfLax.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 11/04/2019.
//

#ifndef SeismicNorm_HopfLax_hpp
#define SeismicNorm_HopfLax_hpp


/** ********* Hopf Lax operator for seismic norms *********

Input : 
- neigh : vertices of the simplex on which to do the optimization.
- l : values of the arrival time at these vertices

Output : 
- (Scalar value) val : estimated arrival time at the center of the stencil.
- (Vector value) w   : coefficients of the geodesic flow.

By construction, w >= 0, and 
sum_i w_i neigh_i is a unit vector, for the norm of interest, the geodesic flow.
*/

// ------ Single vertex --------

template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<1> const& neigh, Vec<1> const& l)
const -> std::pair<ComponentType,Vec<1> > {
	const ScalarType norm0 = Norm(neigh[0]);
	return {norm0+l[0], {1./norm0}};
}



// ---------- Two vertices (edge case) ----------

template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<2> const& neigh, Vec<2> const& l)
const -> std::pair<ComponentType,Vec<2> > {
	static_assert(Dimension>=2);
	// Minimizing a convex function on the interval [0,1]
	VectorType const n01 = neigh[1]-neigh[0];
	ComponentType const l01 = l[1]-l[0];
	assert(l01==l01);
	
	// Check the first vertex
	VectorType const g0 = GradNorm(neigh[0]);
	ComponentType const d0 = g0.ScalarProduct(n01)+l01;
	if(d0>=0) {
		const ScalarType norm0 = g0.ScalarProduct(neigh[0]);
		return {norm0+l[0],Vec<2>{1./norm0,0.}}; }
	//neigh[0].ScalarProduct(g0)
	
	// Check the second vertex
	VectorType const g1 = GradNorm(neigh[1]);
	ComponentType const d1 = g1.ScalarProduct(n01)+l01;
	if(d1<=0) {
		const ScalarType norm1 = g1.ScalarProduct(neigh[1]);
		return {norm1+l[1],Vec<2>{0.,1./norm1}}; }
	//neigh[1].ScalarProduct(g1)
	
	// The minimum is attained in the interior.
	// Make a guess, and call unconstrained optimization
	ComponentType const r = (-d0)/(d1-d0);
	
	if constexpr(Dimension==2) {
		VectorType const neigh_r = neigh[0]*(1-r) + neigh[1]*r;
		ComponentType const l_r = l[0]*(1-r)+l[1]*r;
		return _HopfLax_Face(neigh, l, l_r + Norm(neigh_r));
	} else {
		VectorType g = (1-r)*g0 + r*g1;
		return _HopfLax_Edge(neigh,l,  g);
	}
	
}


template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<3> const& neigh, Vec<3> const& l)
const -> std::pair<ComponentType,Vec<3> > {
	static_assert(Dimension==3);
	
	// We are minimizing a convex function on a triangle
	VectorType const
	n01 = neigh[1]-neigh[0],
	n12 = neigh[2]-neigh[1],
	n20 = neigh[0]-neigh[2];
	
	ComponentType const
	l01 = l[1]-l[0],
	l12 = l[2]-l[1],
	l20 = l[0]-l[2];
	assert(l01==l01 && l12==l12 && l20==l20);
	
	std::array<VectorType,3> const g{
		GradNorm(neigh[0]),
		GradNorm(neigh[1]),
		GradNorm(neigh[2])
	};
	
	ComponentType const
	d01 = g[0].ScalarProduct(n01) + l01,
	d10 = -(g[1].ScalarProduct(n01) + l01),
	d12 = g[1].ScalarProduct(n12) + l12,
	d21 = -(g[2].ScalarProduct(n12) + l12),
	d20 = g[2].ScalarProduct(n20) + l20,
	d02 = -(g[0].ScalarProduct(n20) + l20);
	
	// Is the infimum attained at a vertex ?
	if(d01>=0 && d02>=0){
		const ScalarType norm0 = g[0].ScalarProduct(neigh[0]);
		return {norm0+l[0],VectorType{1./norm0,0.,0.}};}
	
	if(d12>=0 && d10>=0){
		const ScalarType norm1 = g[1].ScalarProduct(neigh[1]);
		return {norm1+l[1],VectorType{0.,1./norm1,0.}};}
	
	if(d20>=0 && d21>=0){
		const ScalarType norm2 = g[2].ScalarProduct(neigh[2]);
		return {norm2+l[2],VectorType{0.,0.,1./norm2}};}
	
	/*		// We use the coordinate chart mapping neigh[0],neigh[1],neigh[2]
	 // onto (0,0), (1,0), (0,1)
	 */
	
	// In the following, ex compute the extrema on each edge, along with the gradient at that position
	
	VectorType edgePos;
	std::array<VectorType,3> edgeGrad;
	std::pair<ComponentType,Vec<3> > result;
	
	
	auto edgeExtremize =
	[&](int i,ComponentType dij, ComponentType dji) -> bool {
		int const j=(i+1)%3, k=(i+2)%3;
		if(dij<0 && dji<0){
			const ComponentType r = dij/(dij+dji);
			VectorType p = (1-r)*g[i]+r*g[j];
			auto const [val,w] =
			//			_HopfLax_Edge({neigh[i],neigh[j]},{l[i],l[j]}); // No guess, seems unreliable
			_HopfLax_Edge({neigh[i],neigh[j]},{l[i],l[j]}, p );
			
			// Alternatively (better guess, but more costly)
			// _HopfLax_Edge({neigh[i],neigh[j]},{l[i],l[j]}, GradNorm((1-r)*neigh[i]+r*neigh[j]) );
			
			assert(abs(p.ScalarProduct(neigh[j]-neigh[i])+l[j]-l[i])<1e-6);
			ComponentType const c =
			p.ScalarProduct(neigh[k]-neigh[j])+l[k]-l[j];
			if(c>=0){
				result.first = val;
				result.second[i]=w[0];
				result.second[j]=w[1];
				result.second[k]=ComponentType(0);
				return true;
			}
			edgePos[i]=w[1]/w.Sum();
			edgeGrad[i]=p;
		} else if(dij>=0){
			edgePos[i]=0;
			edgeGrad[i]=g[i];
		} else { assert(dji>=0);
			edgePos[i]=1;
			edgeGrad[i]=g[j];
		}
		return false;
	};
	
	if(edgeExtremize(0,d01,d10)) return result;
	if(edgeExtremize(1,d12,d21)) return result;
	if(edgeExtremize(2,d20,d02)) return result;
	
	// Now we are building the guess for the face update
	// Find a relation between the gradients
	
	VectorType const a{
		edgeGrad[0].ScalarProduct(n01)+l01,
		edgeGrad[1].ScalarProduct(n01)+l01,
		edgeGrad[2].ScalarProduct(n01)+l01
	}, b{
		edgeGrad[0].ScalarProduct(n12)+l12,
		edgeGrad[1].ScalarProduct(n12)+l12,
		edgeGrad[2].ScalarProduct(n12)+l12
	};
	
	VectorType c = Cross(a,b);
	
	// c should be non-negative, but this is hard to control
	// due to floating point errors in degenerate cases.
	//assert(c.IsNonNegative());
	
	VectorType c2 = VectorType::Constant(0);
	for(int i=0; i<3; ++i){
		int const j=(i+1)%3;
		c2[i]+=c[i]*(1-edgePos[i]);
		c2[j]+=c[i]*edgePos[i];
	}
	
	assert(c2.Sum()!=0);
	c2/=c2.Sum();
	
	// c2 should be non-negative, up to floating point errors
	assert(*std::min_element(c2.begin(),c2.end()) >= -1e-6);
	
	/*
	 if(*std::min_element(c.begin(),c.end()) < -1e-6){
	 std::cout << "In hopf lax "
	 ExportArrayArrow(neigh)
	 ExportVarArrow(l)
	 ExportVarArrow(edgePos)
	 ExportVarArrow(c2)
	 << std::endl;
	 }*/
	
	// Express in terms of the triangle vertices
	
	VectorType const guessPt =
	c2[0]*neigh[0]+c2[1]*neigh[1]+c2[2]*neigh[2];
	
	ComponentType const guessVal = Norm(guessPt)+c2.ScalarProduct(l);
	return _HopfLax_Face(neigh, l, guessVal);
	
}



// Operator associated with interior facet
// Has been shown to fail for anell = 0.020...
template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
_HopfLax_Face(Mat<Dimension> const& neighbors,
			  VectorType const& l, ComponentType const& guess)
const -> std::pair<ComponentType, VectorType> {
	// Create the matrix whose rows are neighbors and invert it
	typedef typename SymmetricMatrixType::MatrixType MatrixType;
	MatrixType V = MatrixType::FromRows(neighbors);
	V = V.Inverse();
	
	typedef Pol1<ComponentType, 1> AffType;
	Vector<AffType,Dimension> z;
	{
		const VectorType u = VectorType::Constant(1);
		const VectorType c1 = V * u;
		const VectorType c0 = V *(guess*u - l);
		for (int i = 0; i < Dimension; ++i) {
			z[i].c[0] = c0[i];
			z[i].c[1] = c1[i];
		}
	}
	
	typedef Pol1<ComponentType, 2*Dimension> DetType;
	
	const DetType det = Constraint<AffType, DetType>(z);
	
	// Apply Newton's root finding algorithm to the determinant,
	// starting from max(L). Yields lambda  (arrival time)
	
	const size_t newtonIt = 8;
	const ComponentType lambda = det.Newton(0., newtonIt);
	// TODO : compare with meaningful value
	assert(std::abs(det(lambda))<1e-6);
	
	typedef DifferentiationType<ComponentType, VectorType> DiffType;
	Vector<DiffType, Dimension> zopt;
	VectorType zopt_raw;
	
	// Find the optimal point in the simplex determined by the neighbors
	for (size_t i = 0; i < Dimension; ++i) {
		zopt_raw[i] = z[i](lambda);
		zopt[i] = DiffType(zopt_raw[i],i);
		// Differentiates along the i-th dimension
	}
	
	DiffType const constraint = Constraint(zopt);
	assert(std::abs(constraint.s) < 1e-6);
	// TODO compare with meaningful value
	
	VectorType weights = V.Transpose()*constraint.v /
		zopt_raw.ScalarProduct(constraint.v);
	// Get position of optimal neighbor
	// (barycentric coefficients, normalized to reproduce unit flow)
//	std::cout << zopt_raw.ScalarProduct(constraint.v) << std::endl;
	
	// TODO : check normalization, ...
	
	// Check that inf is indeed attained in interior, as expected
	// Issue : we can actually get some very small negative values
	
	// Check consistency
	assert((V.Inverse()*GradNorm(V.Transpose().Inverse()*weights)
			+l-(lambda+guess)*VectorType::Constant(1)
			).Norm()<1e-6);
	
	for(int i=0; i<Dimension; ++i){
		assert(weights[i]>=-1e-6);
		weights[i] = std::max(ComponentType(0),weights[i]);
	}
	return { lambda+guess, weights };
}

// Hopf-lax operator associated with edges, in three dimensions
// Rather close to the norm computation
template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
_HopfLax_Edge(Mat<2> const& v, Vec<2> const& l,
			  VectorType & p) // Gradient at optimum. A good guess must be provided
const -> std::pair<ComponentType,Vec<2> > {
	static_assert(Dimension==3, "Edges only arise in dimension three");
	
	using std::sqrt;
	using Sym2 = SymmetricMatrix<ComponentType,2>;
	using Mat32 = Matrix<ComponentType, 3, 2>;
	
	Mat32 const V = Mat32::FromColumns(v);
	
	VectorType const v01 = v[1]-v[0];
	ComponentType const l01 = l[1]-l[0];
	
	const int niter = 12;
	AD2SymType constraint_tmp = AD2Constraint_tmp();
	for (int i = 0; i < niter; ++i) {
		/*
		 // Optimization opportunity : replace autodiff with explicity computations
		 using AD2VectorType = Vector<AD2Type, Dimension>;
		 AD2VectorType pAD;
		 for (int j = 0; j < Dimension; ++j) {
		 pAD[j] = AD2Type(p[j], j);
		 }
		 const AD2Type c = Constraint(pAD);
		 */
		
		const AD2Type c =
		(useMetilOptimization && Dimension==3) ?
		AD2Constraint_Metil(p) :
		AD2Constraint(p,constraint_tmp);
		
		// Finding a line along which the KKT relations lie.
		SymmetricMatrixType const d = c.m.Inverse(); // No minus sign this time
		VectorType const dv01 = d*v01;
		ComponentType const
		c01 = dv01.ScalarProduct(c.v) - (l01+p.ScalarProduct(v01)),
		s0 = v[0].ScalarProduct(dv01),
		s1 = v[1].ScalarProduct(dv01);
		
		Vec<2> pos{s0,s1}; pos *= c01/pos.SquaredNorm();
		Vec<2> const dir{s1,-s0};
		
		VectorType const vpos = V*pos, vdir = V*dir;
		
		// Solving a0+2t*a1+t^2*a2
		ComponentType const
		a0 = d.SquaredNorm(vpos)-d.SquaredNorm(c.v)+2*c.x,
		a1 = d.ScalarProduct(vpos,vdir),
		a2 = d.SquaredNorm(vdir);
		
		ComponentType const delta = a1*a1-a0*a2;
		// Should be positive ... up to floating point error.
		//		assert(delta>=0);
		using std::max;
		
		ComponentType const sdelta = sqrt(max(ComponentType(0.),delta));
		
		/*
		 // Incorrect previous code : choice of quadratic root often wrong
		 ComponentType const orient =
		 d.ScalarProduct(v[0], vdir)/a2;
		 ComponentType const r = (-a1+(orient>0?1:-1)*sdelta)/a2; // sign before sdelta ??
		 
		 Vec<2> const kkt = pos+r*dir;
		 VectorType const w = V*kkt;
		 VectorType const h = d*(w-c.v);
		 
		 p += h;
		 */
		
		ComponentType const r1 = (-a1-sdelta)/a2, r2 = (-a1+sdelta)/a2;
		Vec<2> const kkt1 = pos+r1*dir, kkt2 = pos+r2*dir;
		VectorType const w1 = V*kkt1, w2 = V*kkt2;
		VectorType const h1 = d*(w1 - c.v), h2 = d*(w2-c.v);
		
		const auto b = h1.SquaredNorm() < h2.SquaredNorm();
		using std::abs;
		if constexpr(isSimple) {
			p += b ? h1 : h2;
			// Optimization : adaptive stopping criterion
			if(abs(c.x)<constraintBound_GradNorm) break;
		} else {
#ifdef XSIMD_HPP
			for(int i=0; i<Dimension; ++i) {p[i] = xsimd::select(b,h1[i],h2[i]);}
			if(xsimd::all(abs(c.x)<constraintBound_GradNorm)) break;
#else
			static_assert(isSimple);
#endif
		}
	}
	// Found optimal p.
	ComponentType const val = v[0].ScalarProduct(p)+l[0];
	// Now get the barycentric coefficients
	
	typedef DifferentiationType<ComponentType, VectorType> DiffType;
	Vector<DiffType, Dimension> zopt;
	
	// Find the optimal point in the simplex determined by the neighbors
	for (size_t i = 0; i < Dimension; ++i) {
		zopt[i] = DiffType(p[i],i);
		// Differentiates along the i-th dimension
	}
	
	DiffType const constraint = Constraint(zopt);
	if constexpr(isSimple) {assert(std::abs(constraint.s) < 1e-6);}
	
	// constraint.v should be spanned by the neighbors
	if constexpr(isSimple) {assert(abs(Determinant(v[0], v[1], constraint.v)) < 1e-6);}
	
	const ScalarType normalization = p.ScalarProduct(constraint.v);
	auto const gram = Sym2::EuclideanGram(v);
	Vec<2> weights = // -
	gram.Inverse()*Vec<2>{v[0].ScalarProduct(constraint.v),
		v[1].ScalarProduct(constraint.v)} / normalization;
	if(weights.Sum()<0) weights*=-1;
	
	// Weights should be non-negative .. up to numerical precision
	if constexpr(isSimple) {
		assert(weights[0]/weights.Sum() > -1e-6
			   && weights[1]/weights.Sum() > -1e-6);}
	
	// Should be zero ... up to numerical precision
	if constexpr(isSimple) {assert((GradNorm(V*weights)-p).Norm()<1e-6);}
	
	for(int i=0; i<2; ++i) {weights[i]=std::max(ComponentType(0),weights[i]);}
	return {val,weights};
	
}

#endif /* SeismicNorm_HopfLax_hpp */
