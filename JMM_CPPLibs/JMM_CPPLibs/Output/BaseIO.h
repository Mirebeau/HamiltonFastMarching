// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0


#ifndef BaseIO_h
#define BaseIO_h

#include <map>
#include "IO.h"
#include "JMM_CPPLibs/DataStructures/FromComponentsIterator.h"

/**
 A common base for Python and Mathematica IO.
 */

struct BaseIO : TraitsIO {
    typedef double ScalarType;

    bool HasField(KeyCRef) const;
    bool EraseField(KeyCRef);
    std::string GetString(KeyCRef) const;
    void SetString(KeyCRef, const std::string &);

    BaseIO(){};
    void UsageReport();
protected:
	template<typename TT,typename TC> using FCI=FromComponentsIterator<TT,TC>;
    template<typename T> std::pair<std::vector<DiscreteType>,FCI<T,const ScalarType> > GetDimsPtr(KeyCRef) const;
    template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
    
    struct RawElement;
    std::map<KeyType, RawElement> rawElems;
    
    const RawElement & GetRaw(KeyCRef,bool=true) const;
    RawElement & CreateElement(KeyCRef);
    void SetUsed(KeyCRef);
    SetterTag GetSetter(KeyCRef) const;
    static void StaticSendMsg(bool warn, const std::string & msg) {
        if(warn) std::cout << "***** Warning ! *****\n" << msg << "********************\n";
        else std::cout << msg;}
    virtual void SendMsg(bool warn, const std::string & msg) const {return StaticSendMsg(warn,msg);}
    template<bool,typename> friend struct _Msg;
};

#include "BaseIO.hxx"

#endif /* BaseIO_h */
