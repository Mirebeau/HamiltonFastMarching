// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef ShallowMap_h
#define ShallowMap_h

#include <vector>
#include <cassert>
/*
 This container is made for the (rare) use cases where:
 - one wants to store a sparse array of (usually POD) values.
 - one does not want to pay the storage cost of the full actual array, but an array of indices is fine.
 - one does not want to pay the time complexity cost of an std::map.

 Caution: erase(index) does Not immediately call the destructor of the value stored at this index.
 */

template<typename TIndex, typename TValue, bool hasTrivialConstructor=true> struct ShallowMap {
    typedef TIndex IndexType;
    typedef TValue ValueType;
    
    size_t size() const {return indices.size();}
    void resize(size_t n) {
        for(size_t i=n; i<size(); ++i) if(indices[i]!=BadIndex()) unallocated.push_back(indices[i]);
        indices.resize(n,BadIndex());
        const size_t m = (n*sizeof(IndexType))/(sizeof(ValueType)+sizeof(IndexType));
        reserve_values(m);
    }
    void reserve_values(size_t n){values.reserve(n);unallocated.reserve(n);}
    bool erase(IndexType i) {
        assert(0<=i && i<size());
        IndexType & j = indices[i];
        if(j==BadIndex()) return false;
        unallocated.push_back(j);
            j=BadIndex();
            return true;
    }
    bool find(IndexType i) const {assert(0<=i && i<size()); return indices[i]!=BadIndex();}
    bool empty() const {return values.size()==unallocated.size();}
    const ValueType & operator[](IndexType i) const {
        assert(0<=i && i<size());
        const IndexType j=indices[i];
        return j!=BadIndex() ? values[j] : dummyValue;
    }
    
    ValueType & operator[](IndexType i){
        assert(0<=i && i<size());
        IndexType & j=indices[i];
        if(j!=BadIndex()) return values[j];
        else {
            if(unallocated.empty()){
                j=(IndexType)values.size();
                Construct::Emplace(values);
//                values.emplace_back();
            } else {
                j=unallocated.back();
                unallocated.pop_back();
                values[j]=ValueType();
            }
            return values[j];
        }
    }
    template<typename ...Args>
    ValueType & emplace(IndexType i, Args&&...args){
        IndexType &j=indices[i];
        if(j==BadIndex()){
            if(unallocated.empty()){
                j=values.size();
                values.emplace_back(std::move(args...));
                return values[j];
            } else {
                j=unallocated.back();
                unallocated.pop_back();
            }
        }
        values[j]=ValueType(std::move(args...));
        return values[j];
    }
protected:
    std::vector<IndexType> indices;
    std::vector<IndexType> unallocated;
    
    std::vector<ValueType> values;
    ValueType dummyValue;
    static constexpr IndexType BadIndex() {return std::numeric_limits<IndexType>::max();}
    template<bool hasTrivialConstructor_,typename Dummy> struct Construct_;
	typedef Construct_<hasTrivialConstructor, void> Construct;
};

template<typename TI, typename TV, bool hTC> template<typename Dummy>
struct ShallowMap<TI,TV,hTC>::Construct_<true,Dummy>{
    static void Emplace(std::vector<TV> & values){values.emplace_back();}
};

template<typename TI, typename TV, bool hTC> template<typename Dummy>
struct ShallowMap<TI,TV,hTC>::Construct_<false,Dummy>{
    static void Emplace(std::vector<TV> & values){assert(false);}
};



#endif /* ShallowMap_h */
