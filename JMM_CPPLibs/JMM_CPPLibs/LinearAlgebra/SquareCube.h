// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_SquareCube_h
#define AmongConvex2_SquareCube_h


template<int i, typename T>
inline T pow(T t){
    if(i<0) return T(1)/pow<-i,T>(t);
    if(i==0) return T(1);
    if(i==1) return t;
    const int j=i%2;
    const T s = pow< (i>>1),T>(t);
    if(j==0) return s*s;
    else return t*s*s;
}


template<typename T>
inline T square(T t){return t*=t;}

template<typename T>
inline T cube(T t){return t*=square(t);}

// non-negative modulo
template<typename T> inline T PosMod(T a, T n){
    assert(n>0);
    const T b = a%n;
    return b>=0 ? b : b+n;
}

template<typename T>
inline T fPosMod(T a, T n){
    assert(n>0);
    const T b = fmod(a,n);
    return b>=0 ? b : b+n;
}

#endif
