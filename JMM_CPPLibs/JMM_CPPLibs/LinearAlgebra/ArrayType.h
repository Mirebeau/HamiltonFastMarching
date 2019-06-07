// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef ArrayType_h
#define ArrayType_h

/*
 A very basic array type.
 Row major ordering. (Changing from Column major on 20/11/2017)
 */

#include <vector>
#include <algorithm>
#include "PointType.h"
#include "../Macros/RedeclareTypes.h"


namespace LinearAlgebra {

template<typename TComponent, size_t VDim, typename TDiscrete=int>
struct Array : public std::vector<TComponent> {
    typedef TComponent ComponentType;
    typedef TDiscrete DiscreteType;
    static const size_t Dimension = VDim;
    typedef Point<DiscreteType, Dimension> IndexType;
    
    typedef std::vector<ComponentType> Superclass;
    Array(Superclass && super_, IndexType dims_):Superclass(std::move(super_)),dims(dims_){};
    Array(const Superclass & super_,IndexType dims_):Superclass(super_),dims(dims_){};
    Array(){};
    
    IndexType dims=IndexType::Constant(0);
    bool CheckDims() const {return dims.Product()==this->size();}
    bool InRange(const IndexType & index) const {
        for(size_t i=0; i<Dimension; ++i) if(index[i]<0 || index[i]>=dims[i]) return false;
        return true;}
	Redeclare2Types(Superclass, reference, const_reference);
    reference operator() (const IndexType & index) {assert(CheckDims() && InRange(index));
        return this->operator[](Convert(index));}
    const_reference operator() (const IndexType & index) const {assert(CheckDims() && InRange(index));
        return this->operator[](Convert(index));}
    
    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (std::ostream & os, const Array & a){a.PrintSelf(os); return os;}

    DiscreteType Convert(const IndexType & index) const { // Row Major
        if(Dimension==0) return 0;
        assert(InRange(index));
        DiscreteType result = index[0];
        for(int i=1; i<Dimension; ++i)
            result = result*dims[i]+index[i];
        return result;
    }
    
    IndexType Convert(DiscreteType n) const { // Row major
        if(Dimension==0) return IndexType();
        assert(0<=n && n<dims.Product());
        IndexType result;
        for(int i=int(Dimension)-1; i>0; --i){
            result[i]=n%dims[i];
            n/=dims[i];
        }
        result[0]=n;
        return result;
    }
    
    template<typename F> auto
    Transform(const F & f) const -> Array<decltype(f(this->operator[](0))),VDim,TDiscrete>  {
        typedef decltype(f(this->operator[](0))) TC;
        Array<TC,VDim,TDiscrete> result;
        result.dims = dims;
        result.resize(this->size());
        std::transform(this->begin(),this->end(),result.begin(),f);
        return result;
    }
    
    template<typename TC> Array<TC,VDim,TDiscrete> Cast() const {
        auto caster = [](const TComponent & a)->TC{return TC(a);};
        return this->Transform(caster);
    }
};
    
template<typename TC, size_t VD, typename TD> void
Array<TC,VD,TD>::PrintSelf(std::ostream & os) const {
    if(!CheckDims()) {
        os << "{\"Dimensions " << dims << "inconsistent with size " << this->size() << "\"}";
        return;
    }
    IndexType mods;
    mods[0]=dims[0];
    for(int i=1; i<Dimension;++i)
        mods[i]=dims[i]*mods[i-1];
    
    for(DiscreteType n=0; n<this->size();){
        for(int i=0; i<Dimension; ++i){
            if(n%mods[i]==0) os << "{";
            else break;
        }
        os << this->operator[](n);
        ++n;
        for(int i=0; i<Dimension; ++i){
            if(n%mods[i]==0) os << "}";
            else break;
        }
        if(n<this->size()) os << ",";
        else break;
    }
}

    
}

#endif /* ArrayType_h */
