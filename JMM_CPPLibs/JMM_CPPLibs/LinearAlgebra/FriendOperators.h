// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_FriendOperators_h
#define AmongConvex2_FriendOperators_h

// Operators for a vector space V on a field K. (or a module on a ring, but be careful with divisions.)
// K elements are passed by value, since we assume that it is a machine type (e.g. double, int)

namespace LinearAlgebra {
    
template<class K> class multiplicative {
    friend K operator * (K k, const K & l){k*=l; return k;}
    friend K operator / (K k, const K & l){k/=l; return k;}
};

// Define u+=v, u-=v, u*=k, u/=k
template <class V, class K> class vector_space {
    friend V operator + (V u, const V & v) {u+=v; return u;}
    friend V operator - (V u, const V & v) {u-=v; return u;}
    friend V operator * (V u, const K & k)  {u*=k; return u;}
    friend V operator * (const K & k, V u)  {u*=k; return u;}
    friend V operator / (V u, const K & k)  {u/=k; return u;}
};

// define p+=v, p-=v
template <class P, class V> class offsettable {
    friend P operator +(P p, const V & v){p+=v; return p;}
    friend P operator +(const V & v, P p){p+=v; return p;}
    friend P operator -(P p, const V & v){p-=v; return p;}
};

/*
// Natural implementation, given below, causes some compilers (VS15) to allocate space
// in order to disambiguate the addresses of the two empty classes.
// The C++20 attribute [[no_unique_address]] had no effect when tested

template <class P, class V, class K> class affine_space :
offsettable< P, V>,
vector_space<V, K>
{};


// We thus switch to a more verbose implementation below
*/

template<class P, class V, class K> class affine_space :
vector_space<V, K> {
	// copied from offsettable< P, V>
	friend P operator +(P p, const V & v) { p += v; return p; }
	friend P operator +(const V & v, P p) { p += v; return p; }
	friend P operator -(P p, const V & v) { p -= v; return p; }
};

// define a<b
template <class A> class totally_ordered {
    friend bool operator > (const A &a, const A &b){return (b<a);};
    friend bool operator <= (const A &a, const A &b){return !(b<a);}
    friend bool operator >= (const A &a, const A &b){return !(a<b);}
};

// define a<b and a>b
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
    
}

#endif
