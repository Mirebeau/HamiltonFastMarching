// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef VirtualConstIterator_h
#define VirtualConstIterator_h

template<typename TElement>
struct VirtualConstIterator {
    typedef TElement ElementType;
    virtual ~VirtualConstIterator(){};
    
    virtual void operator++() = 0;
    bool IsAtEnd() const {return status==AtEnd;}
    const ElementType & operator*() const {return current;}
    VirtualConstIterator(int _status, const ElementType & _current):status(_status),current(_current){}
protected:
    enum {AtEnd=0};
    int status;
    ElementType current;
};

#endif /* VirtualConstIterator_h */
