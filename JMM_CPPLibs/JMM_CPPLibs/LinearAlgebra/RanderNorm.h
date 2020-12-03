// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef RanderNorm_LinearAlgebra_h
#define RanderNorm_LinearAlgebra_h

#include "VectorType.h"
#include "SymmetricMatrixType.h"

namespace LinearAlgebra {
    
template<typename TScalar, size_t VDimension>
struct RanderNorm {
    
    using ScalarType = TScalar;
    static const size_t Dimension = VDimension;

    using VectorType = Vector<ScalarType,Dimension>;
    using SymmetricMatrixType = SymmetricMatrix<ScalarType,Dimension>;
    
    SymmetricMatrixType m;
    VectorType w;
    
    ScalarType Norm(const VectorType & v) const {return m.Norm(v) - w.ScalarProduct(v);}
    const VectorType Gradient(const VectorType & v) const {return m*v/m.Norm(v) - w;}
	
    const RanderNorm DualNorm() const;
    bool IsDefinite() const {return m.Inverse().SquaredNorm(w)<1;}
};
    
template<typename TS, size_t VD>
const RanderNorm<TS,VD> RanderNorm<TS,VD>::DualNorm() const {
    const SymmetricMatrixType s = (m-SymmetricMatrixType::RankOneTensor(w)).Inverse();
    const VectorType omega = s*w;
    return RanderNorm{s*(1+w.ScalarProduct(omega)),omega};
}

template<typename TS, size_t VD>
std::ostream & operator << (std::ostream & os, const RanderNorm<TS,VD> & norm){
    return os << "{" << norm.m << "," << norm.w << "}";
}
	
} // namespace Linear Algebra

#endif
