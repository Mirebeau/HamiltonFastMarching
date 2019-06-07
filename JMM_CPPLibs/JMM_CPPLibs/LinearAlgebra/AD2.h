#ifndef AD2_h
#define AD2_h

/*
 This file implements second order forward automatic differentiation, in arbitrary but fixed dimension.
 */

#include <cmath>
#include "SymmetricMatrixType.h"

namespace LinearAlgebra {

template<typename TAD, typename TComponent> struct AD2Operators;
	
// second order automatic differentiation class
template<typename TComponent, size_t VDimension>
struct AD2 : AD2Operators<AD2<TComponent,VDimension>,TComponent> {
	typedef TComponent ComponentType;
	static const size_t Dimension = VDimension;
	typedef LinearAlgebra::SymmetricMatrix<ComponentType, Dimension> SymmetricMatrixType;
	typedef LinearAlgebra::Vector<ComponentType, Dimension> VectorType;

	ComponentType x;
	VectorType v;
	SymmetricMatrixType m;

	// Represents the Taylor expansion : x + eps.v + 1/2*eps.m.eps

	AD2() {};

	AD2(ComponentType x_, int i) {
		assert(0 <= i && i < Dimension);
		x = x_;
		v = VectorType::Constant(0);
		v[i] = 1;
		m = SymmetricMatrixType::Zero();
	}

	AD2(ComponentType x_) {
		x = x_;
		v = VectorType::Constant(0);
		m = SymmetricMatrixType::Zero();
	}

	void operator+= (AD2 const& a) {
		x += a.x;
		v += a.v;
		m += a.m;
	};

	void operator+= (ComponentType const& a) {
		x += a;
	};

	void operator-= (AD2 const& a) {
		x -= a.x;
		v -= a.v;
		m -= a.m;
	};

	void operator-= (ComponentType const& a) {
		x -= a;
	};

	AD2 operator - () const {
		AD2 a; 
		a.x = -x;
		a.v = -v;
		a.m = -m;
		return a;
	}

	void operator*= (AD2 const& a) {
		int const d = Dimension;
		m = x*(a.m)+(a.x)*m;
		// On effectue m+= v*(a.v)+(a.v)*v
		for (int i = 0; i < d; ++i) {
			for (int j = i; j < d; ++j) {
				m(i, j) += v[i] * a.v[j] + a.v[i] * v[j];
			}
		}
		v = v*(a.x) + x*(a.v);
		x *= a.x;
	};

	void operator*= (ComponentType const& a) {
		x *= a;
		v *= a;
		m *= a;
	};


};

template<typename TC,size_t VD>
std::ostream & operator << (std::ostream & os, AD2<TC, VD> const & a) {
	os << "{" << a.x << "," << a.v << "," << a.m << "}";
	return os;
}

template<typename TC,size_t VD>
AD2<TC, VD> sqrt(AD2<TC, VD> a) {
	assert(a.x >= 0);
	using std::sqrt;
	TC y = sqrt(a.x);
	a.m = a.m / (TC(2.) * y);
	for (int i = 0; i < VD; ++i) {
		for (int j = 0; j <= i; ++j) {
			a.m(i, j) -= a.v[i] * a.v[j] / (TC(4.) * a.x*y);
		}
	}
	a.v = TC(1. / 2.) * a.v / sqrt(a.x);
	a.x = y;
	return a;
}


/*
// Derivative of a two argument function
template<typename funct>
AD2<double,2> deriv(funct const & f, double x, double y) {
	typedef AD2<double, 2> ADType;
	return f(ADType(x, 0), ADType(y, 1));
}
*/


#include "Implementation/AD2.hpp"
	
} // namespace LineaAlgebra

#endif // AD2_h
