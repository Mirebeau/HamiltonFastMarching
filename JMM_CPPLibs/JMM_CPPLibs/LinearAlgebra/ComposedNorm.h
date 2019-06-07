//
//  NormComposition.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 27/03/2019.
//

#ifndef ComposedNorm_h
#define ComposedNorm_h
namespace LinearAlgebra {
	
	/// Implements the composition of a norm with a (linear) transformation
template <typename TNorm, typename TTrans, size_t TDim = TNorm::Dimension>
struct ComposedNorm {
	static const int Dimension = TDim;
	using NormType = TNorm;
	using TransformType = TTrans;
	using VectorType = typename NormType::VectorType;
	using ComponentType = typename NormType::ComponentType;
	NormType m;
	TransformType a;
	ComponentType Norm(VectorType x) const { return m.Norm(a*x); }
	VectorType Gradient(VectorType x) const { return a.Transpose()*m.Gradient(a*x); }
	
	template<size_t n> using Vec = typename NormType::template Vec<n>; //LinearAlgebra::Vector<ComponentType,n>;
	template<size_t n> using Mat = typename NormType::template Mat<n>; //std::array<VectorType,n>;
	
	// Hopf-Lax update operator
	template<size_t n> std::pair<ComponentType, Vec<n> >
	HopfLax(Mat<n> const & neigh, Vec<n> const & l) const {
		return m.HopfLax(Transform(neigh),l);}

	// Optimized variants, avoiding multiple recomputations. Additional arguments are caching data.
	std::pair<ComponentType, Vec<1> > HopfLax(Mat<1> const & neigh, Vec<1> const & l,
											  VectorType & cacheOut) const {
		return m.HopfLax(Transform(neigh),l,cacheOut);}
	std::pair<ComponentType, Vec<2> > HopfLax(Mat<2> const & neigh, Vec<2> const & l,
											  Mat<2> const & cacheIn, VectorType & cacheOut) const {
		return m.HopfLax(Transform(neigh),l,cacheIn,cacheOut);}
	std::pair<ComponentType, Vec<3> > HopfLax(Mat<3> const & neigh, Vec<3> const & l,
		Mat<3> const & vertexCacheIn, Mat<3> const & edgeCacheIn, Vec<3> const & edgeCacheIn2) const {
		return m.HopfLax(Transform(neigh),l,vertexCacheIn,edgeCacheIn,edgeCacheIn2);}
	
protected:
	template<size_t n> Mat<n> Transform(Mat<n> const & neigh) const {
		Mat<n> result;
		for(int i=0; i<n; ++i) result[i] = a * neigh[i];
		return result;}
};
	
	
/// Implements a special kind of linear transformation, determined by a single vector
template<typename TVec, bool VTrans = false>
struct TopographicTransform {
	using VectorType = TVec;
	static const bool transposed = VTrans;
	VectorType v;
	VectorType operator *(VectorType const & x) const {
		VectorType result = x;
		if (!transposed) {result.back() += v.ScalarProduct(x);}
		else {result += v*x.back();}
		return result;
	}
	
	using TransposedTransformType = TopographicTransform<VectorType,!transposed>;
	TransposedTransformType Transpose() const {
		return TransposedTransformType{ v };}
};
	
}


#endif /* NormComposition_h */
