#ifndef Pol1_h
#define Pol1_h

/*
 This file implements polynomials of an arbitrary but fixed degree.
 */

#include "VectorType.h"
#include "AD2.h"

namespace LinearAlgebra {
	
template<typename TPol1, typename TComponent> struct Pol1Operators;

// classe pour les polynomes, avec opérateurs définis par héritage
template<typename TComponent, size_t VDegree>
struct Pol1 :
Pol1Operators<Pol1<TComponent, VDegree>, TComponent>
{
	typedef TComponent ComponentType;
	static const size_t Degree = VDegree;
	typedef Vector<ComponentType, Degree+1> DataType;

	Pol1() = default;
	Pol1(ComponentType const & x) {
		c = DataType::Constant(0);
		c[0] = x;
	}
	explicit Pol1(DataType const & c_):c(c_){};

	DataType c; //Polynomial coefficients
	
	void operator += (Pol1 const& q) {
		c += q.c;
	}
	void operator += (ComponentType const& x) {
		c[0] += x;
	}
	void operator -= (Pol1 const& q) {
		c -= q.c;
	}
	void operator -= (ComponentType const& x) {
		c[0] -= x;
	}
	void operator *= (ComponentType const& x) {
		c *= x;
	}
	void operator /= (ComponentType const& x) {
		assert(x!=0);
		c/=x;
	}
	Pol1 operator-() const {
		Pol1 neg;
		neg.c=-c;
		return neg;
	}
	template<size_t Degree2> 
	Pol1<ComponentType, Degree + Degree2> operator * (Pol1<ComponentType, Degree2> const & q) const {
		Pol1<ComponentType, Degree + Degree2> r; // result
		std::fill(r.c.begin(), r.c.end(), ComponentType(0));
		for (size_t i = 0; i <= Degree; ++i) {
			for (size_t j = 0; j <= Degree2; ++j) {
				r.c[i + j] += c[i] * q.c[j];
			}
		}
		return r;
	}

	ComponentType operator() (ComponentType const& x) const {
		ComponentType xi = 1, result = 0;
		for (size_t i = 0; i <= Degree; ++i) {
			result += c[i] * xi;
			xi *= x;
		}
		return result;
	}

	Pol1<ComponentType, Degree - 1> Deriv() const {
		Pol1<ComponentType, Degree - 1> q;
		for (size_t i = 1; i <= Degree; ++i) {
			q.c[i - 1] = i * c[i];
		}
		return q;
	}
	
	// Apply Newton iterations so as to find a root, starting from given guess.
	ComponentType Newton(ComponentType x, size_t it) const {
		const Pol1 & self = *this;
		const auto deriv = self.Deriv();
		ComponentType result = x;
		for (size_t i = 0; i < it; ++i) {
			assert(deriv(result) != 0);
			result -= self(result) / deriv(result);
		}
/*		std::cout << "resultat de la methode de Newton (x, P(x)) : " << result << ", " << self(result) << std::endl;
		std::cout << "polynome et depart : " << self << ", " << x << std::endl;*/
		return result;
	}
	
	/*
	// Apply quadratic Newton-like iterations,
	// so as to find a root larger than guess
	ComponentType Newton2(ComponentType x, size_t it) const {
		const Pol1 & self = *this;
		const auto deriv = self.Deriv();
		const auto deriv2 = deriv.Deriv();
		ComponentType r=x; // returned value
		for(size_t i=0; i<it; ++i){
			ComponentType const a = deriv2(r)/2, b=deriv(r)/2, c=self(r);
			ComponentType const delta = b*b-a*c;
			assert(delta>=0);
			using std::sqrt;
			ComponentType const sdelta = sqrt(delta);
			assert(a!=0);
			using std::abs;
			ComponentType const r0=r-b/a, h=sdelta/abs(a);
			r = r0-h;
			if(r<x) {r = r0+h;}
			std::cout
			ExportVarArrow(a)
			ExportVarArrow(b)
			ExportVarArrow(c)
			ExportVarArrow(r0-h)
			ExportVarArrow(r0+h)
			ExportVarArrow(r)
			<< std::endl;
		}
		return r;
	}
*/
};

template<typename TC, size_t VD>
std::ostream & operator << (std::ostream & os, Pol1<TC, VD> const & pol) {
	os << pol.c;
	return os;
}

template<typename TPol1, typename TComponent>
struct Pol1Operators {
	typedef TPol1 Pol1;
	typedef TComponent Comp;
	friend Pol1 operator + (Pol1 a, Pol1 const & b) { a += b; return a; }
	friend Pol1 operator + (Pol1 a, Comp const & b) { a += b; return a; }
	friend Pol1 operator + (Comp const & a, Pol1 b) { b += a; return b; }
	
	friend Pol1 operator - (Pol1 a, Pol1 const & b) { a -= b; return a; }
	friend Pol1 operator - (Pol1 a, Comp const & b) { a -= b; return a; }
	friend Pol1 operator - (Comp const & a, Pol1 b) { b -= a; return -b; }
	
	friend Pol1 operator * (Pol1 a, Comp const & b) { a *= b; return a; }
	friend Pol1 operator * (Comp const & a, Pol1 b) { b *= a; return b; }
	
	friend Pol1 operator / (Pol1 a, Comp const & b) { a/=b; return a;}
};

}

#endif // Pol1_h
