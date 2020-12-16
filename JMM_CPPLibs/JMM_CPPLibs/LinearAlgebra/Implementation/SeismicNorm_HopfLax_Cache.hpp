//
//  SeismicNorm_HopfLax_Cache.hpp
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 11/04/2019.
//

#ifndef SeismicNorm_HopfLax_Cache_hpp
#define SeismicNorm_HopfLax_Cache_hpp


// ********* Hopf Lax operator for seismic norms *********

/*
 The following code is extremely close to the non-caching operators, and thus
 violates the DRY rule.
 TODO : express the original operators in terms of the caching ones, or suppress them.
 */

// ------ Single vertex --------


template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<1> const& neigh, Vec<1> const& l,
		VectorType & cache) // Cached data : Gradient of norm at argument
const -> std::pair<ComponentType,Vec<1> > {
	const VectorType grad0 = GradNorm(neigh[0]);
	cache = grad0; // Caching for reuse
	const ComponentType norm0 = grad0.ScalarProduct(neigh[0]);
	return {norm0+l[0], {1.}};
}


// --- Two vertices (edge case) ----

template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<2> const& neigh, Vec<2> const& l,
		Mat<2> const& cacheIn, // Gradients at the vertices
		VectorType& cacheOut) // Gradient at the optimal point (returned in dimension 3 only)
const -> std::pair<ComponentType,Vec<2> > {
	static_assert(Dimension>=2);
	
	// Minimizing a convex function on the interval [0,1]
	VectorType const n01 = neigh[1]-neigh[0];
	ComponentType const l01 = l[1]-l[0];
	assert(l01==l01);
	
	// Check the first vertex
	VectorType const g0 = cacheIn[0]; //GradNorm(neigh[0]);
	assert((g0-GradNorm(neigh[0])).Norm()<1e-6*g0.Norm());
	ComponentType const d0 = g0.ScalarProduct(n01)+l01;
	if(d0>=0) {cacheOut = g0; return {g0.ScalarProduct(neigh[0])+l[0],Vec<2>{1.,0.}}; }
	
	// Check the second vertex
	VectorType const g1 = cacheIn[1]; //GradNorm(neigh[1]);
	assert((g1-GradNorm(neigh[1])).Norm()<1e-6*g1.Norm());
	ComponentType const d1 = g1.ScalarProduct(n01)+l01;
	if(d1<=0) {cacheOut = g1; return {g1.ScalarProduct(neigh[1])+l[1],Vec<2>{0.,1.}}; }
	
	// The minimum is attained in the interior.
	// Make a guess, and call unconstrained optimization
	ComponentType const r = (-d0)/(d1-d0);
	
	if constexpr(Dimension==2) {
		// cacheOut untouched by design
		VectorType const neigh_r = neigh[0]*(1-r) + neigh[1]*r;
		ComponentType const l_r = l[0]*(1-r)+l[1]*r;
		return _HopfLax_Face(neigh, l, l_r + Norm(neigh_r));
	} else {
		// cacheOut is gradient
		cacheOut = (1-r)*g0+r*g1;
		return _HopfLax_Edge(neigh,l, cacheOut);
	}
}

// --- Three vertices ---


template<typename TC, size_t VD> auto SeismicNorm<TC, VD>::
HopfLax(Mat<3> const & neigh, Vec<3> const& l,
		Mat<3> const & vertexCacheIn, // Gradients at the vertices
		Mat<3> const & edgeCacheIn, // Gradients at the optimal position on the edge
		Vec<3> const & edgeCacheIn2 // Position of optimum on the edge
		)
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
	
	/* // Previous code with recomputation
	std::array<VectorType,3> const g{
		GradNorm(neigh[0]),
		GradNorm(neigh[1]),
		GradNorm(neigh[2])
	};
	*/
	const auto & g = vertexCacheIn;
	
	ComponentType const
	d01 = g[0].ScalarProduct(n01) + l01,
	d10 = -(g[1].ScalarProduct(n01) + l01),
	d12 = g[1].ScalarProduct(n12) + l12,
	d21 = -(g[2].ScalarProduct(n12) + l12),
	d20 = g[2].ScalarProduct(n20) + l20,
	d02 = -(g[0].ScalarProduct(n20) + l20);
	
	// Is the infimum attained at a vertex ?
	if(d01>=0 && d02>=0){
		return {g[0].ScalarProduct(neigh[0])+l[0],VectorType{1.,0.,0.}};}
	
	if(d12>=0 && d10>=0){
		return {g[1].ScalarProduct(neigh[1])+l[1],VectorType{0.,1.,0.}};}
	
	if(d20>=0 && d21>=0){
		return {g[2].ScalarProduct(neigh[2])+l[2],VectorType{0.,0.,1.}};}
	
	/*		// We use the coordinate chart mapping neigh[0],neigh[1],neigh[2]
	 // onto (0,0), (1,0), (0,1)
	 */
	
	// In the following, ex compute the extrema on each edge, along with the gradient at that position

	/* // Previous code with recomputation
	 VectorType edgePos;
	 std::array<VectorType,3> edgeGrad;
	 */

	/* The i-th edge has vertices neigh[i] and neigh[j], j=i+1
	 Minimum along this edge is attained at
	 (1-edgePos[i])*neigh[i]+edgePos[i]*neigh[j],
	 The gradient of the norm at this minimum is edgeGrad[i]
	 */
	const auto & edgeGrad = edgeCacheIn;
	const auto & edgePos = edgeCacheIn2;
	std::pair<ComponentType,Vec<3> > result;
	
	
	auto edgeExtremize =
	[&](int i,ComponentType dij, ComponentType dji) -> bool {
		int const j=(i+1)%3, k=(i+2)%3;
		if(dij<0 && dji<0){ // Min of edge is interior to edge
	
			// Previous code with recomputation
//			const ComponentType r = dij/(dij+dji);
//			auto const [val,w,p] =
//			_HopfLax_Edge({neigh[i],neigh[j]},{l[i],l[j]}, (1-r)*g[i]+r*g[j] );
			
			// Alternatively (better guess, but more costly)
			// _HopfLax_Edge({neigh[i],neigh[j]},{l[i],l[j]}, GradNorm((1-r)*neigh[i]+r*neigh[j]) );
			
			const VectorType & p = edgeGrad[i];
			
			assert(abs(p.ScalarProduct(neigh[j]-neigh[i])+l[j]-l[i])<1e-6);
			ComponentType const c =
			p.ScalarProduct(neigh[k]-neigh[j])+l[k]-l[j];
			if(c>=0){ // Min of edge is global min
				const ComponentType wj = edgePos[i];
				const ComponentType wi = 1-wj;
				
/*				// Previous code with recomputation
 				result.first = val;
				result.second[i]=w[0];
				result.second[j]=w[1];
				result.second[k]=ComponentType(0);*/
				
				result.first = (wi*neigh[i]+wj*neigh[j]).ScalarProduct(p) + wi*l[i]+wj*l[j];
				result.second[i]=wi;
				result.second[j]=wj;
				result.second[k]=ComponentType(0);
				
				return true;
			}
			// Old code with recompute
//			edgePos[i]=w[1]/w.Sum();
//			edgeGrad[i]=p;
		} else if(dij>=0){
			// Old code with recompute
//			edgePos[i]=0;
//			edgeGrad[i]=g[i];
			assert(edgePos[i]==0);
			assert(edgeGrad[i]==g[i]);
		} else { assert(dji>=0);
			// Old code with recompute
//			edgePos[i]=1;
//			edgeGrad[i]=g[j];
			assert(edgePos[i]==1);
			
/*			std::cout << " In HL3 "
			ExportVarArrow(edgeGrad[i])
			ExportVarArrow(g[j])
			ExportArrayArrow(g)
			ExportArrayArrow(edgeGrad)
			ExportArrayArrow(edgePos)
			ExportVarArrow(i)
			<< std::endl;*/
			 
			assert(edgeGrad[i]==g[j]);
		}
		return false;
	};
	
	if(edgeExtremize(0,d01,d10)) return result;
	if(edgeExtremize(1,d12,d21)) return result;
	if(edgeExtremize(2,d20,d02)) return result;
	
	// Now we are building the guess for the face update
	// c2 will store the barycentric coordinates of the guess w.r.t neighbors
	VectorType c2 = VectorType::Constant(0);
	
	// First, check wether one of the vertices is a good guess
	const ScalarType tol=1e-3;
	if     (edgePos[1]>=1-tol && edgePos[2]<=tol){c2[2]=1;} 
	else if(edgePos[2]>=1-tol && edgePos[0]<=tol){c2[0]=1;}
	else if(edgePos[0]>=1-tol && edgePos[1]<=tol){c2[1]=1;}
	else {
	/* Otherwise, determine the guess by finding a relation between the gradients
	computed at the edge minimizers.
	If the minimized function was quadratic, this would be the minimum.*/

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

	for(int i=0; i<3; ++i){
		int const j=(i+1)%3;
		c2[i]+=c[i]*(1-edgePos[i]);
		c2[j]+=c[i]*edgePos[i];
	}
	
	assert(c2.Sum()!=0);
	c2/=c2.Sum();
	
	// c2 should be non-negative, up to floating point errors
	
	/*
	 if(*std::min_element(c2.begin(),c2.end()) < -1e-6){
		 std::cout << "In hopf lax "
		 ExportArrayArrow(neigh)
		 ExportVarArrow(l)
		 ExportVarArrow(edgePos)
		 ExportVarArrow(c2)
		 ExportVarArrow(c)
		 ExportVarArrow(a)
		 ExportVarArrow(b)
		 ExportArrayArrow(edgeGrad)
		 ExportVarArrow(l01)
		 ExportVarArrow(l12) << "\n"
		 ExportVarArrow(d01)
		 ExportVarArrow(d10)
		 ExportVarArrow(d02)
		 ExportVarArrow(d20)
		 ExportVarArrow(d12)
		 ExportVarArrow(d21)
		 << std::endl;
	 }*/

	assert(*std::min_element(c2.begin(),c2.end()) >= -1e-6);

	// Express in terms of the triangle vertices
	}
		
	const VectorType guessPt = c2[0]*neigh[0]+c2[1]*neigh[1]+c2[2]*neigh[2];
	ComponentType const guessVal = Norm(guessPt)+c2.ScalarProduct(l);
	return _HopfLax_Face(neigh, l, guessVal);
	
}


#endif /* SeismicNorm_HopfLax_Cache_hpp */


