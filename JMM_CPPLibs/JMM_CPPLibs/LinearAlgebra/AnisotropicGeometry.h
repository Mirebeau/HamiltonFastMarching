//
//  AnisotropicGeometry.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 05/03/2019.
//

#ifndef AnisotropicGeometry_h
#define AnisotropicGeometry_h

namespace LinearAlgebra {

// A generalized cosine, symmetric in its arguments
template<typename NormType, typename VectorType,
typename ComponentType = typename VectorType::ComponentType>
ComponentType CosAngle(const NormType & norm, const VectorType & u, const VectorType & v){
	const VectorType gu = norm.Gradient(u), gv = norm.Gradient(v);
	ComponentType const nu = gu.ScalarProduct(u), nv=gv.ScalarProduct(v);
	assert(nu>0 && nv>0);
	using std::min;
	return min(gu.ScalarProduct(v)/nv, gv.ScalarProduct(u)/nu);
}

template<typename NormType, typename VectorType>
bool IsAcute(const NormType & norm, const VectorType & u, const VectorType & v){
	const VectorType gu = norm.Gradient(u), gv = norm.Gradient(v);
	return gu.ScalarProduct(v)>=0 && gv.ScalarProduct(u)>=0;
}

}
#endif /* AnisotropicGeometry_h */
