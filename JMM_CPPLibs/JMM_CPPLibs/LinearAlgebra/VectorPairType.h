// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef VectorPairType_h
#define VectorPairType_h
#include "FriendOperators.h"
#include "../DataStructures/GetComponent.h"

namespace LinearAlgebra {
template<typename TV0, typename TV1, typename TComponent = typename TV0::ComponentType>
struct VectorPair : 
vector_space< VectorPair<TV0, TV1, TComponent>, TComponent>,
std::pair<TV0,TV1>
 {
    typedef std::pair<TV0,TV1> Superclass;
    using Superclass::Superclass;
    typedef TComponent ComponentType;
    
    // Define u+=v, u-=v, u*=k, u/=k
    VectorPair & operator+=(const VectorPair & u){this->first+=u.first; this->second+=u.second; return *this;}
    VectorPair & operator-=(const VectorPair & u){this->first-=u.first; this->second-=u.second; return *this;}
    VectorPair & operator*=(const ComponentType & a){this->first*=a; this->second*=a; return *this;}
    VectorPair & operator/=(const ComponentType & a){this->first/=a; this->second/=a; return *this;}
    
    VectorPair operator -() const {return VectorPair(-this->first, -this->second);}
};

}

template<typename TA, typename TB, typename C>
struct GetComponent<LinearAlgebra::VectorPair<TA,TB>, C>
: GetComponent<std::pair<TA,TB>, C> {};
#endif /* VectorPairType_h */
