// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_FriendOperators_h
#define AmongConvex2_FriendOperators_h


namespace LinearAlgebra {

// --------- Types holding the result of an operation -------------

template<typename S,typename T> using product_t = decltype(S()*T());
template<typename S,typename T> using sum_t = decltype(S()+T());
template<typename S,typename T> using algebra_t =
std::enable_if_t<std::is_same_v<product_t<S,T>, sum_t<S,T> >, sum_t<S,T> >;

// ---------- Operator overloading ----------

/* The following structures are used to defined some c++ operators in terms of others.
Operators for a vector space V on a field K.
(or a module on a ring, but be careful with divisions.)*/

// Requires *=, /=
template<class K> class multiplicative {
    friend K operator * (K k, const K & l){k*=l; return k;}
    friend K operator / (K k, const K & l){k/=l; return k;}
};

// Requires u+=v, u-=v, u*=k, u/=k
template <class V, class K> class vector_space {
    friend V operator + (V u, const V & v) {u+=v; return u;}
    friend V operator - (V u, const V & v) {u-=v; return u;}
    friend V operator * (V u, const K & k)  {u*=k; return u;}
    friend V operator * (const K & k, V u)  {u*=k; return u;}
    friend V operator / (V u, const K & k)  {u/=k; return u;}
};

// Requires p+=v, p-=v
template <class P, class V> class offsettable {
    friend P operator +(P p, const V & v){p+=v; return p;}
    friend P operator +(const V & v, P p){p+=v; return p;}
    friend P operator -(P p, const V & v){p-=v; return p;}
};


// Requires a<b
template <class A> class totally_ordered {
    friend bool operator >  (const A &a, const A &b){return (b<a);};
    friend bool operator <= (const A &a, const A &b){return !(b<a);}
    friend bool operator >= (const A &a, const A &b){return !(a<b);}
};

// Requires a<b and a>b
template <class A,class B> class totally_ordered2 {
	friend bool operator <= (const A & a, const B & b){return !(a>b);}
	friend bool operator >= (const A & a, const B & b){return !(a<b);}
	
	friend bool operator <  (const B & b, const A & a){return a>b;}
	friend bool operator >  (const B & b, const A & a){return a<b;}
	friend bool operator <= (const B & b, const A & a){return a>=b;}
	friend bool operator >= (const B & b, const A & a){return a<=b;}
};
// define a==b
template <class A> class equality_comparable {
    friend bool operator != (const A &a, const A &b){return !(a==b);}
};

// ---------------- Multiple operators ---------------

/* The most compact, and safe, implementation of multiple additional operators
 would be based on multiple inheritance. However, this causes some compilers (VS15) to
 allocate space in order to disambiguate the addresses of the multiple superclasses.
 We thus switch to another more verbose implementation.
 The C++20 attribute [[no_unique_address]] had no effect when tested*/


/* Equivalent to the following multiple inheritance.
template <class P, class V, class K> class affine_space :
offsettable< P, V>, vector_space<V, K> {}; */
template<class P, class V, class K> class affine_space :
vector_space<V, K> {
	// copied from offsettable< P, V>
	friend P operator +(P p, const V & v) { p += v; return p; }
	friend P operator +(const V & v, P p) { p += v; return p; }
	friend P operator -(P p, const V & v) { p -= v; return p; }
};


/* Equivalent to the following multiple inheritance.
 
vector_space< DifferentiationType<TScalar,TVector>, TScalar>,
offsettable<DifferentiationType<TScalar,TVector>, TScalar>,
totally_ordered<DifferentiationType<TScalar, TVector> >,
totally_ordered2<DifferentiationType<TScalar, TVector>, TScalar>,
equality_comparable<DifferentiationType<TScalar, TVector> >,
multiplicative<DifferentiationType<TScalar, TVector> > */
template<typename TDifferentiation, typename TScalar> struct DifferentiationTypeOperators {
	using AD = TDifferentiation;
	using K = TScalar;
	friend AD operator + (AD a, AD const & b) { a += b; return a; }
	friend AD operator + (AD a, K const & b)  { a += b; return a; }
	friend AD operator + (K const & a, AD b)  { b += a; return b; }

	friend AD operator - (AD a, AD const & b) { a -= b; return a; }
	friend AD operator - (AD a, K const & b)  { a -= b; return a; }
	friend AD operator - (K const & a, AD b)  { b -= a; return -b; }

	friend AD operator * (AD a, AD const & b) { a *= b; return a; }
	friend AD operator * (AD a, K const & b)  { a *= b; return a; }
	friend AD operator * (K const & a, AD b)  { b *= a; return b; }
	
	friend AD operator / (const AD & a, const AD & b){return a*b.Inverse();}
	friend AD operator / (const K & a,  const AD & b){return a*b.Inverse();}
	friend AD operator / (AD a, const K & b) {a/=b; return a;}
	friend void operator /= (AD & a, const AD & b){a*=b.Inverse();}

	// Totally ordered
	friend bool operator >  (AD const & a, AD const & b) { return b < a; }
	friend bool operator >= (AD const & a, AD const & b) { return !(a < b); }
	friend bool operator <= (AD const & a, AD const & b) { return !(a > b); }

	// Totally ordered 2
	friend bool operator <= (const AD & a, const K & b) { return !(a > b); }
	friend bool operator >= (const AD & a, const K & b) { return !(a < b); }

	friend bool operator <  (const K & b, const AD & a) { return a > b; }
	friend bool operator >  (const K & b, const AD & a) { return a < b; }
	friend bool operator <= (const K & b, const AD & a) { return a >= b; }
	friend bool operator >= (const K & b, const AD & a) { return a <= b; }

	friend bool operator != (const AD &a, const AD &b) { return !(a == b); }
};

}

#endif
