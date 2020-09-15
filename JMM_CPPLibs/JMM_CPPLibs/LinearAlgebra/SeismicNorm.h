#ifndef SeismicNorm_h
#define SeismicNorm_h

/* This file implements norms arising seismology, defined in terms of the fourths order Hooke elasticity tensor of the domain.
 
 // TODO : measure of anellipticity
 */

#include "SymmetricMatrixType.h"
#include "DifferentiationType.h"
#include "AD2.h"
#include "Pol1.h"
#include "AnisotropicGeometry.h"
#include "ScalarTraits.h"

namespace LinearAlgebra {

template<typename TComponent, size_t VDimension>
struct SeismicNorm {
	using ComponentType = TComponent;
	static const int Dimension = VDimension;
	using SymmetricMatrixType =  SymmetricMatrix<ComponentType, Dimension>;
	using VectorType = Vector<ComponentType, Dimension>;
	static const int VoigtDimension = (Dimension*(Dimension + 1)) / 2;
	using HookeTensorType = SymmetricMatrix<ComponentType, VoigtDimension>;
	HookeTensorType hookeTensor;

	SeismicNorm() = default;
	SeismicNorm(SymmetricMatrixType const &); 	/// Constructor for elliptical norms
	SeismicNorm(HookeTensorType const & ht):hookeTensor(ht){};

	ComponentType const& Coefficient(int i, int j, int k, int l) const {
		return hookeTensor(VoigtIndex(i, j), VoigtIndex(k, l));
	};

	ComponentType& Coefficient(int i, int j, int k, int l) {
		return hookeTensor(VoigtIndex(i, j), VoigtIndex(k, l));
	};

	int VoigtIndex(int, int) const;
	std::pair<int, int> VoigtIndexInv(int) const;

	// Tensor contraction c_{ijkl} u_j u_l
	template<typename T, typename TOut=T>
	SymmetricMatrix<TOut, Dimension>
	Contract(Vector<T, Dimension> const &) const;
	// Returns det(lambda*Id - c_{ijkl} u_j u_l)
	template<typename T, typename TOut=T, typename S=ComponentType> TOut
	Constraint(Vector<T, Dimension> const &, S const & = S(1.)) const;
	
	// Primal norm, defining the metric
	VectorType GradNorm(VectorType const &) const;
	ComponentType Norm(VectorType const& q) const {
		return q.IsNull() ? 0. : GradNorm(q).ScalarProduct(q);
	}
	VectorType Gradient(VectorType const & q) const {return GradNorm(q);}
	
	// Dual norm, appearing in the eikonal equation
	VectorType GradDualNorm(VectorType const &) const;
	template<typename T> T
	DualNorm(Vector<T,Dimension> const &) const;

	// Hopf-Lax update operator, involved in eikonal solvers
	template<size_t n> using Vec = LinearAlgebra::Vector<ComponentType,n>;
	template<size_t n> using Mat = std::array<VectorType,n>;
	std::pair<ComponentType, Vec<1> > HopfLax(Mat<1> const&, Vec<1> const &) const;
	std::pair<ComponentType, Vec<2> > HopfLax(Mat<2> const&, Vec<2> const &) const;
	std::pair<ComponentType, Vec<3> > HopfLax(Mat<3> const&, Vec<3> const &) const;
	
	// Optimized variants, avoiding multiple recomputations.
	// Additional arguments are caching data (gradient at the optimal point in normalized coords)
	std::pair<ComponentType, Vec<1> > HopfLax(Mat<1> const&, Vec<1> const &,
											  VectorType &) const;
	std::pair<ComponentType, Vec<2> > HopfLax(Mat<2> const&, Vec<2> const &,
											  Mat<2> const &, VectorType &) const;
	std::pair<ComponentType, Vec<3> > HopfLax(Mat<3> const&, Vec<3> const &,
											  Mat<3> const &, Mat<3> const &, Vec<3> const &) const;

	// Geometric measures
	ComponentType Anellipticity() const; // Low value diagnoses non-smoothness
	SymmetricMatrixType EllipticApproximation() const;
protected:
	// "Unsafe" HopfLax operators, assuming the assuming minimizer is interior,
	// and requiring a good guess
	std::pair<ComponentType, VectorType>
	_HopfLax_Face(Mat<Dimension> const&, VectorType const&, ComponentType const&) const;
	
	// Three dimensional only variant, devoted to edges
	std::pair<ComponentType, Vec<2> > _HopfLax_Edge(Mat<2> const&, Vec<2> const&, VectorType &) const;
	
	constexpr static const int nIterMax_GradNorm = 8;
	
	constexpr static const bool isSimple = is_simple_v<ComponentType>;
#ifdef XSIMD_HPP
	using ScalarType = xsimd::revert_simd_type<ComponentType>;
#else
	using ScalarType = ComponentType;
#endif
	constexpr static const ScalarType
	constraintBound_GradNorm = 1000*std::numeric_limits<ScalarType>::epsilon(),
	constraintBound_HopfLax_Edge = 1000*std::numeric_limits<ScalarType>::epsilon();
	/*
	constexpr static const ComponentType
	constraintBound_GradNorm = 1000*std::numeric_limits<ComponentType>::epsilon(),
	constraintBound_HopfLax_Edge = 1000*std::numeric_limits<ComponentType>::epsilon();
	*/
	// Internal methods for accelerating GradNorm and _HopfLax_Edge computations
	using AD2Type = AD2<ComponentType, Dimension>;
	using AD2SymType = LinearAlgebra::SymmetricMatrix<AD2Type, Dimension>;
	AD2SymType AD2Constraint_tmp() const;
	AD2Type AD2Constraint(VectorType const &, AD2SymType &) const;
	AD2Type AD2Constraint_zero(AD2SymType const &) const;
	AD2Type AD2Constraint_Metil(VectorType const &) const;
	static const bool useMetilOptimization = true;
};


template <typename TC,size_t VD>
std::ostream & operator << (std::ostream & os, const SeismicNorm<TC, VD> & norm) {
	os << norm.hookeTensor; 
	return os;
}

	
#include "Implementation/SeismicNorm.hpp"
#include "Implementation/SeismicNorm_HopfLax.hpp"
#include "Implementation/SeismicNorm_HopfLax_Cache.hpp"
#include "Implementation/SeismicNorm_MetilConstraint.hpp"


} // namespace LinearAlgebra

#endif // SeismicNorm_h 
