//
//  TTINorm.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 02/01/2020.
//

#ifndef TTINorm_h
#define TTINorm_h

#include "SymmetricMatrixType.h"
#include "JMM_CPPLibs/Macros/DependentFalse.h"
#include "AD2.h"
#include "BasisReduction.h"
#include "SquareCube.h"

namespace LinearAlgebra {

/**
 Tilted transversly anisotropic norm.
 A class of norms used in seismic traveltime tomography.
 Less general than the SeismicNorm type, which allows for a generic Hooke tensor.
 Allows for a causal Eulerian discretization.
 
 The dual unit ball is defined by
 < linear,p > + (1/2)< p,quadratic,p > = 1
 where p is the vector containing the squares of transform*p0.
 
 */
template<typename TScalar, int VDimension>
struct TTINorm {
	using ScalarType = TScalar;
	static const int Dimension = VDimension;
	
	using LinearType = Vector<ScalarType,2>;
	using QuadraticType = SymmetricMatrix<ScalarType, 2>;
	
	using VectorType = Vector<ScalarType,Dimension>;
	using TransformType = Matrix<ScalarType, Dimension, Dimension>;
	using SymmetricMatrixType = SymmetricMatrix<ScalarType, Dimension>;
	LinearType linear;
	QuadraticType quadratic;
	TransformType transform;
	
	/// Level set for the dual unit ball.
	template<typename T> T Level(const Vector<T,Dimension> &) const;
	VectorType Gradient(const VectorType &) const;
	ScalarType Norm(const VectorType &) const;
	
	/// Gauss-Siedel update operator
	static const int SymDimension = SymmetricMatrixType::InternalDimension;
	using NeighborValuesType = std::array<ScalarType,SymDimension>;
	
	/// Properties related to the expression of the TTI norm as an extremum of Riemannian norms.
	struct Properties {
		int optimDirection; /// wether a sup or an inf.
		ScalarType tMin, tMax; /// Diagonal anisotropy is (1-t,t) or (1-t,t,t)
		PrintSelfMacro(Properties);
	};
	Properties Props() const;
	ScalarType smallest_positive_root(ScalarType, ScalarType) const;
	
	using DiscreteType = int;
	using ReductionType =
	LinearAlgebra::BasisReduction<ScalarType, DiscreteType, Dimension>;
	Redeclare3Types(ReductionType,DiscreteVectorType,SuperbaseType,SellingPath)
	SellingPath Selling(ScalarType) const;
	template<typename T> T Multiplier(const T &) const;
	template<typename T> T UpdateValue(const T &, const NeighborValuesType &,
									   const SellingPath &) const;
	PrintSelfMacro(TTINorm);
};

#include "Implementation/TTINorm.hxx"
}
#endif /* TTINorm_h */
