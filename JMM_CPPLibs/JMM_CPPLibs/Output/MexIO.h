// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MexIO_h
#define MexIO_h

#include "mex.h"
#include <sstream>
#include <chrono>

#define __IgnoreTrailingSingletonDimensions 1
#include "BaseIO.h"


    
struct MexIO : TraitsIO {
    typedef double ScalarType;
    
    template<bool warn> struct Msg_ {
        std::ostringstream oss;
        template<typename T> Msg_ & operator << (const T & t){oss << t; return (*this);}
        ~Msg_();
    };
    
    bool HasField(KeyCRef) const;
    const mxArray * HasMxField(KeyCRef) const;
    std::string GetString(KeyCRef) const;
    void SetString(KeyCRef, std::string);
    
    MexIO(const mxArray *, mxArray **);
    MexIO(const BaseIO &) = delete;
    static std::chrono::time_point<std::chrono::system_clock> top;
    void UsageReport();
protected:
    template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;
    template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
    mutable std::set<KeyType> unused, defaulted, exported;
    void SetField(KeyCRef, mxArray *);
    
    const mxArray * mxInput;
    mxArray ** mxOutput;
    static void StaticSendMsg(bool warn, std::string msg);
    void SendMsg(bool warn, const std::string & msg) const {return StaticSendMsg(warn,msg);}
    template<bool,typename> friend struct _Msg;
};
std::chrono::time_point<std::chrono::system_clock> MexIO::top = std::chrono::system_clock::now();
#include "MexIO.hxx"



#endif /* MexIO_h */
