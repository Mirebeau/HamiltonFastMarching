// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Licensed under the Apache License, Version 2.0, WITHOUT ANY WARRANTY, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef CappedVector_h
#define CappedVector_h

#include <array>
#include <iostream>

template<typename T, size_t N>
struct CappedVector :  std::array<T,N> {
    typedef std::array<T,N> Superclass;

    using Superclass::begin;
    using Superclass::cbegin;
    
    using typename Superclass::value_type;
    using typename Superclass::size_type;
    using typename Superclass::reference;
    using typename Superclass::const_reference;
    using typename Superclass::iterator;
    using typename Superclass::const_iterator;

    iterator end() {return _end;}
    const_iterator end() const {return _end;}
    const_iterator cend() const {return _end;}
    
    iterator rbegin() {return Superclass::rbegin() + std::distance(_end,Superclass::end());}
    const_iterator rbegin() const {return Superclass::rbegin() + std::distance(_end,Superclass::end());}
    const_iterator crbegin() const {return rbegin();}
    
    reference operator[](size_type n){assert(n<=max_size()); return Superclass::operator[](n);}
    const_reference operator[](size_type n) const {assert(n<=max_size()); return Superclass::operator[](n);}
    reference at(size_type n){if(n>max_size()) throw std::logic_error("Capped vector index exceeds size"); return Superclass::at(n);}
    const_reference at(size_type n) const {if(n>max_size()) throw std::logic_error("Capped vector index exceeds size"); return Superclass::at(n);}
    
    reference front(){assert(!empty()); return Superclass::front();}
    const_reference front() const {assert(!empty()); return Superclass::front();}
    reference back(){assert(!empty()); auto __end = _end; return *(--__end);}
    const_reference back() const {assert(!empty()); auto __end = _end; return *(--__end);}
    
    void push_back(const_reference value){assert(size()<max_size()); *_end=value; ++_end;}
    void pop_back(){assert(!empty()); --_end;}
    
    bool empty() const {return cbegin()==cend();}
    size_type size() const {return std::distance(cbegin(),cend());}
    static constexpr size_type max_size() {return N;}
    
    void clear(){_end=begin();}
    void resize(size_type n, const_reference val=value_type()){assert(n<=max_size());
        if(n<=size()) _end=begin()+n;
        else {iterator __end = _end; _end=begin()+n; std::fill(__end, _end, val);}}
    void reserve(size_type n){if(n>N) throw std::logic_error("Capped vector reserve cannot exceed max_size");}
    void fill(const_reference value){std::fill(begin(),end(),value);}
    
    CappedVector(){clear();}
    CappedVector(std::initializer_list<T> l){resize(l.size()); std::copy(l.begin(), l.end(), begin());}
    CappedVector(const CappedVector & other){resize(other.size()); std::copy(other.begin(),other.end(),begin());}
    CappedVector & operator=(const CappedVector & other){resize(other.size()); std::copy(other.begin(),other.end(),begin()); return*this;}
protected:
    iterator _end;
};

template<typename T, size_t N>
std::ostream & operator << (std::ostream & os, const CappedVector<T, N> & v){
    os << "{";
    if(!v.empty()){
        auto it = v.begin();
        os << *it; ++it;
        for(;it!=v.end();++it)
            os << "," << *it;
    }
    os << "}";
	return os;
}

#endif /* CappedVector_h */
