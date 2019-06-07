// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef RangeAccessor_h
#define RangeAccessor_h

template<typename TPointer> struct RangeAccessor
: std::pair<TPointer,TPointer> {
    using Pointer = TPointer;
	using Superclass = std::pair<Pointer,Pointer>;
	RangeAccessor(Pointer begin, Pointer end):Superclass{begin,end}{};
    
    Pointer begin() const {return this->first;}
    Pointer end() const {return this->second;}
    size_t  size() const {return end() - begin();}
    
    typedef decltype(*Pointer()) value_type;
    value_type & operator[](size_t i) {assert(i<size()); return begin()[i];}
    const value_type & operator[](size_t i) const {assert(i<size()); return begin()[i];}
    value_type & back() {assert(size()>0); return *(--end());}
    const value_type & back() const {assert(size()>0); return *(--end());}
    value_type & front() {assert(size()>0); return *begin();}
    const value_type & front() const {assert(size()>0); return *begin();}
};

#endif /* RangeAccessor_h */
