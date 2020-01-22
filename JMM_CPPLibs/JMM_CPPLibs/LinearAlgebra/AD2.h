#ifndef AD2_h
#define AD2_h

/*
 This file implements second order forward automatic differentiation,
 in arbitrary but fixed dimension.
 An auxiliary function implements a single step of Sequential-Quadratic-Programming.
 */

#include <cmath>

#include "SymmetricMatrixType.h"
#include "UtilsAD.h"

namespace LinearAlgebra {
	
// second order automatic differentiation class
template<typename TComponent, size_t VDimension> struct AD2 :
DifferentiationTypeOperators<AD2<TComponent,VDimension>,TComponent> {
	static constexpr size_t Dimension = VDimension;
	using ComponentType =  TComponent;
	using VectorType = LinearAlgebra::Vector<ComponentType, Dimension>;
	using SymmetricMatrixType = LinearAlgebra::SymmetricMatrix<ComponentType, Dimension> ;
	using Sym = SymmetricMatrixType;
	
	ComponentType x;
	VectorType v;
	SymmetricMatrixType m;

	// Represents the Taylor expansion : x + eps.v + 1/2*eps.m.eps

	AD2() {};

	AD2(const ComponentType & x_, int i) {
		assert(0 <= i && i < Dimension);
		x = x_;
		v = VectorType::Constant(0);
		v[i] = 1;
		m = SymmetricMatrixType::Zero();
	}

	AD2(const ComponentType & x_) {
		x = x_;
		v = VectorType::Constant(0);
		m = SymmetricMatrixType::Zero();
	}
	
	AD2(const ComponentType & x_, const VectorType & v_, const SymmetricMatrixType & m_)
	:x(x_),v(v_),m(m_){};
	
	// --------- Addition and substraction ----------
	
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

	// ---------- Multiplication and division ----------
	void operator*= (AD2 const& a) {
		m = x*(a.m)+(a.x)*m + Sym::Outer2(v,a.v);
		v = v*(a.x) + x*(a.v);
		x *= a.x;
	};

	void operator*= (ComponentType const& a) {
		x *= a;
		v *= a;
		m *= a;
	};
	
	AD2 Inverse() const {
		assert(x!=0);
		const ComponentType xi = 1./x, xi2=xi*xi;
		return AD2{xi, -v*xi2, xi2*((2*xi)*Sym::OuterSelf(v) - m)};
	}
		
	void operator/= (ComponentType const& a) {
		x /= a;
		v /= a;
		m /= a;
	};
	// ----------- Comparison -----------
	bool operator<(AD2 const& a) const {return x<a.x;}
	bool operator<(ComponentType const& a){return x<a;}
	
	// ------------ Special functions ---------
	friend AD2 sqrt(const AD2 & a){
		assert(a>=0.);
		using std::sqrt;
		const ComponentType x=a.x, xs=sqrt(x), xi=1/x;
		const VectorType w = (0.5*xi)*a.v;
		return xs*AD2{1.,w, 0.5*(xi*a.m - Sym::OuterSelf(w))};
	}

	static Vector<AD2,Dimension> Perturbation(const VectorType & v){
		Vector<AD2,Dimension> delta;
		for(int i=0; i<Dimension;++i){delta[i]=AD2(v[i],i);}
		return delta;
	}
	
	/* ------ Single step of Sequential Quadratic Programming -------
	 Maximize <q,p> subject to constraint >=0, differentiated at p. Usage :
	 for(int i=0; i<nIter; ++i){p+=Level(Perturbation(p)).SQP(q);}
	 */
	VectorType SQP(const VectorType & q) const {
		const AD2 & c = *this;
		const SymmetricMatrixType d = c.m.Inverse();
		const VectorType dv = d*c.v, dq = d*q;
		const ComponentType num = dv.ScalarProduct(c.v) - 2.*c.x;
		const ComponentType den = dq.ScalarProduct(q);
		assert(num*den >= 0);
		using std::sqrt;
		const ComponentType lambda = -sqrt(num / den);
		const VectorType h = lambda*dq - dv;
		return h;
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

template<typename TScalar, size_t VDimension>
struct UtilsAD<AD2<TScalar,VDimension> >{
	using ScalarType = TScalar;
	using ADType = LinearAlgebra::AD2<TScalar,VDimension>;
	static ScalarType Scalar(const ADType & t){return t.x;}
};

/*
// Derivative of a two argument function
template<typename funct>
AD2<double,2> deriv(funct const & f, double x, double y) {
	typedef AD2<double, 2> ADType;
	return f(ADType(x, 0), ADType(y, 1));
}
*/


#include "Implementation/AD2.hpp"
	
} // namespace LinearAlgebra

#endif // AD2_h
