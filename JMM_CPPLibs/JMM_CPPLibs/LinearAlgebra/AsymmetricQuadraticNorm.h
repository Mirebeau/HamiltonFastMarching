// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AsymmetricQuadraticNorm_h
#define AsymmetricQuadraticNorm_h

/*
 
 Implements the norm v.M.v+max(w.v,0)^2
 
 */

#include "SymmetricMatrixType.h"
#include "AnisotropicGeometry.h"

namespace LinearAlgebra {
    
template<typename TScalar, size_t VDimension>
struct AsymmetricQuadraticNorm {
    
    typedef TScalar ScalarType;
    static const size_t Dimension = VDimension;
    
    typedef Vector<ScalarType, Dimension> VectorType;
    typedef SymmetricMatrix<ScalarType, Dimension> SymmetricMatrixType;
    
    SymmetricMatrixType m;
    VectorType w;
    
    ScalarType Norm(const VectorType & v) const {
        const ScalarType scal = ScalNeg(v);
        return sqrt(m.SquaredNorm(v)+scal*scal);
    }
    
    ///Positive multiple of gradient
    VectorType MGrad(const VectorType & v) const {return m*v+ScalNeg(v)*w;}
    VectorType Gradient(const VectorType & v) const {
        const VectorType g = MGrad(v);
        const ScalarType squaredNorm = g.ScalarProduct(v);
        return g/sqrt(squaredNorm);
    }
	
    const AsymmetricQuadraticNorm DualNorm() const {
		assert(IsDefinite());
		const SymmetricMatrixType M = (m+SymmetricMatrixType::RankOneTensor(w)).Inverse();
		const VectorType wInv = m.CGSolve(w);
		using std::sqrt;
		const VectorType W = -wInv/sqrt(1+w.ScalarProduct(wInv));
		return AsymmetricQuadraticNorm{M,W};
	}
    bool IsDefinite() const {return m.IsPositive();}
    
protected:
    ScalarType ScalNeg(const VectorType & v) const {return std::min(w.ScalarProduct(v),ScalarType(0));}
};
    
template<typename TS, size_t VD>
std::ostream & operator << (std::ostream & os, const AsymmetricQuadraticNorm<TS,VD> & norm){
    return os << "{" << norm.m << "," << norm.w << "}";
}


template<typename ScalarType, size_t Dimension, typename VectorType>
bool IsAcute(const AsymmetricQuadraticNorm<ScalarType,Dimension> & norm,
			 const VectorType & u, const VectorType & v) {
	const ScalarType
	muv = norm.m.ScalarProduct(u,v), wu = norm.w.ScalarProduct(u), wv=norm.w.ScalarProduct(v);
	using std::min;
	return muv+min(wu*min(wv,0.),wv*min(wu,0.)) >= 0.;
}

} // namespace LineaAlgebra


#endif /* AsymmetricQuadraticNorm_h */
