// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

// This file implements the norm defined as v->\| d * v\|, where d is a fixed vector
// and d * v denotes coefficientwise multiplication.


#ifndef DiagonalNorm_h
#define DiagonalNorm_h

#include "VectorType.h"

namespace LinearAlgebra {

template<typename TScalar, size_t VDimension>
struct DiagonalNorm {
	using ScalarType = TScalar;
	static const size_t Dimension = VDimension;
	
    using VectorType = Vector<ScalarType,Dimension>;
	VectorType d;
	
	ScalarType Norm(const VectorType & v) const {
		ScalarType s=0.;
		for(int i=0; i<Dimension; ++i) s += square(d[i]*v[i]);
		using std::sqrt;
		return sqrt(s);
	}
	const VectorType Gradient(const VectorType & v) const {
		VectorType g;
		for(int i=0; i<Dimension; ++i) g[i] = square(d[i])*v[i];
		using std::sqrt;
		return g/sqrt(g.ScalarProduct(v));
	}
	const DiagonalNorm DualNorm() const {
		DiagonalNorm dual;
		for(int i=0; i<Dimension; ++i) dual.d[i] = 1./d[i];
		return dual;
	}
	bool IsDefinite() const {return d.IsPositive();}
};

template<typename TS, size_t VD>
std::ostream & operator << (std::ostream & os, const DiagonalNorm<TS,VD> & norm){
    return os << norm.d;
}

}


#endif /* DiagonalNorm_h */
