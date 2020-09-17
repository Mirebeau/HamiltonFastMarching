// Copyright 2020 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AsymRanderNorm_h
#define AsymRanderNorm_h

/* This file implements "Asymmetric Rander norms",
 a class of Finslerian, strongly asymmetric norm
 proposed by Da Chen. The norm is determined by a positive definite
 matrix and three vector fields. It is defined as
 N(x) = sqrt( <x,M x> + <u,x>_^2 + <v,x>_^2) - <w,x>
where a_ = min(0,a)
 The HFM fast marching method only makes this norm available in two dimensions,
 due to the difficulty of developing a causal scheme in higher dimensions.
 */

#include "VectorType.h"
#include "SymmetricMatrixType.h"
#include "VectorPairType.h"
#include "HopfLaxMinimize.h"

namespace LinearAlgebra {

template<typename TScalar, size_t VDimension>
struct AsymRanderNorm {
	using ScalarType = TScalar;
	static const size_t Dimension = VDimension;
	
	template<size_t n> using Vec = Vector<ScalarType,n>;
	using VectorType = Vec<Dimension>;
	using PointType = typename VectorType::PointType;
	using SymmetricMatrixType = SymmetricMatrix<ScalarType,Dimension>;
	using AddableSourceType = VectorPair<SymmetricMatrixType,Vec<3*Dimension> >;
	
	SymmetricMatrixType m;
	VectorType u,v,w;
	AsymRanderNorm(){}
	AsymRanderNorm(const AddableSourceType & s):m(s.first){
		auto it = s.second.begin();
		for(int i=0;i<Dimension; ++i){u[i]=*it++;}
		for(int i=0;i<Dimension; ++i){v[i]=*it++;}
		for(int i=0;i<Dimension; ++i){w[i]=*it++;}
	}
		
	ScalarType Norm(const VectorType & x) const {
		using std::sqrt; using std::min;
		const ScalarType
		ux_  = min(0.,u.ScalarProduct(x)),
		vx_  = min(0.,v.ScalarProduct(x)),
		wx   = w.ScalarProduct(x),
		mxx  = m.SquaredNorm(x);
		
		return sqrt(mxx + ux_*ux_ + vx_*vx_) - wx;
	}
	const VectorType Gradient(const VectorType & x) const {
		using std::sqrt; using std::min;
		const VectorType mx = m*x;
		const ScalarType
		ux_  = min(0.,u.ScalarProduct(x)),
		vx_  = min(0.,v.ScalarProduct(x)),
		mxx  = mx.ScalarProduct(x);
		return (mx + u*ux_ + v*vx_)/sqrt(mxx + ux_*ux_ + vx_*vx_) - w;
	}
	
	// DualNorm -> Pruposedly not implemented, since structure is not preserved
	
	/// This is a necessary condition, but not a sufficient condition. The vector w must also be small enough.
	bool IsDefinite(){return m.IsDefinite();}
	
	// Out : center value, geodesic flow coefs. In : neighbors, neighbor values.
	std::pair<ScalarType,VectorType>
	HopfLax(std::array<VectorType,2> const &, VectorType const &) const;
	// Vertex adapted variant, somewhat dummy
	std::pair<ScalarType,Vec<1> >
	HopfLax(std::array<VectorType,1> const &, Vec<1> const &) const;
	
	/// Apply a prescribed diagonal change of coordinates
	void Rescale(const PointType & h){
		for(int i=0; i<Dimension; ++i){
			for(int j=0; j<=i; ++j){
				m(i,j) *= h[i]*h[j];
			}
		}
		for(int i=0; i<Dimension; ++i){
			u[i]*=h[i];
			v[i]*=h[i];
			w[i]*=h[i];
		}
	}
}; // AsymRanderNorm


template<typename TScalar, size_t VDimension> auto
AsymRanderNorm<TScalar,VDimension>::
HopfLax(std::array<VectorType,2> const & neigh, VectorType const & l)
const -> std::pair<ScalarType,VectorType> {
			
	// ------ Optimize in the interior ------
	// Perform change of variables, to reduce to the case of canonical neighbors [1,0],[0,1]
	const SymmetricMatrixType m = this->m.Gram(neigh);
	const VectorType u = {this->u.ScalarProduct(neigh[0]),this->u.ScalarProduct(neigh[1])};
	const VectorType v = {this->v.ScalarProduct(neigh[0]),this->v.ScalarProduct(neigh[1])};
	const VectorType w = {this->w.ScalarProduct(neigh[0]),this->w.ScalarProduct(neigh[1])};
	
	// Introduce modified neighbor values, taking into account the 'drift' w.
	const VectorType lw = l-w;
	const SymmetricMatrixType
	uu = SymmetricMatrixType::OuterSelf(u),
	vv = SymmetricMatrixType::OuterSelf(v);
	
	const ScalarType inf = std::numeric_limits<ScalarType>::infinity();
	
	// Solve with more or less asymmetric quadratic terms.
	for(int eu=0; eu<=1; ++eu){
		for(int ev=0; ev<=1; ++ev){
			SymmetricMatrixType M = m;
			if(eu) M += uu;
			if(ev) M += vv;
			
			VectorType flow;
			const ScalarType val = HopfLaxMinimize(M,lw,flow);
			if(val==inf) continue;
			// Check that the flow has the correct scalar products
			if(eu != (u.ScalarProduct(flow)<=0)) continue;
			if(ev != (v.ScalarProduct(flow)<=0)) continue;
			
			flow /= (1-w.ScalarProduct(flow));
			return {val,flow};
		}
	}
	
	// ------ Check the vertices of the simplex ------
	// Optimization opportunity : compute norms in the transformed coordinate system above
	const PointType norms = {Norm(neigh[0]), Norm(neigh[1])};
	const PointType vals = norms+l;
	if(vals[0]<vals[1]) return {vals[0],{1./norms[0],0.}};
	else {return {vals[1],{0.,1./norms[1]}};}
}

template<typename TScalar, size_t VDimension> auto
AsymRanderNorm<TScalar,VDimension>::
HopfLax(std::array<VectorType,1> const & neigh, Vec<1> const & l)
const -> std::pair<ScalarType,Vec<1> > {
	const ScalarType norm = Norm(neigh[0]);
	return {l[0]+norm,1./norm};
}

	
} // namespace LinearAlgebra
#endif /* AsymRanderNorm_h */
