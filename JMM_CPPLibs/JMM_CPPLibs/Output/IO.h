// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef IO_h
#define IO_h

/** A simple import/export interface. 
 Note that most of the complexity is due to Matlab(R)'s habit of tampering with trailing singleton dimensions, and of permuting the first two coordinate axes in images.
 
 The code WILL throw an exception, with a descriptive error message, whenever the user input is invalid. Make sure that your user input code is exception safe.
 */

#include <string>
#include <set>
#include "../Macros/Exception.h"
#include "../LinearAlgebra/ArrayType.h"
#include "EnumToString.h"
#include "../Macros/RedeclareTypes.h"


#ifndef __IgnoreTrailingSingletonDimensions
#define __IgnoreTrailingSingletonDimensions 0
#endif



struct TraitsIO {
    template<typename TC, size_t VD> using Array = LinearAlgebra::Array<TC,VD>;
    typedef std::string KeyType;
    typedef std::string const & KeyCRef;
    typedef Array<double,1>::DiscreteType DiscreteType;
    TraitsIO(const TraitsIO &) = delete;
    TraitsIO(){};

	enum class SetterTag {User,Compute,Unknown};
    enum class ArrayOrdering {RowMajor, ColumnMajor, YXZ_RowMajor, YXZ_ColumnMajor, Ignore};

	SetterTag currentSetter = SetterTag::User;
    int verbosity=1;
    ArrayOrdering arrayOrdering = ArrayOrdering::RowMajor;
protected:
    template<size_t d> using DimType = LinearAlgebra::Point<DiscreteType,d>;
    mutable std::set<KeyType> unused, defaulted, visitedUnset;
};

template<bool warn, typename IO> struct _Msg {
    std::ostringstream oss;
    const IO * pio = nullptr;
    _Msg(){};
    _Msg(const IO * _pio):pio(_pio){};
    ~_Msg(){
        const std::string & msg = oss.str();
        if(pio) pio->SendMsg(warn,msg);
        else IO::StaticSendMsg(warn,msg);}
    template<typename T> _Msg & operator << (const T & t){oss << t; return *this;}
};

template<> char const * enumStrings<TraitsIO::ArrayOrdering>::data[] = {"RowMajor", "ColumnMajor", "YXZ_RowMajor", "YXZ_ColumnMajor"};
template<> char const * enumStrings<TraitsIO::SetterTag>::data[] = {"User","Compute","Unknown"};

/* ---- INFO : DO NOT DELETE -----
 Base must inherit TraitsIO (or redefine similar traits) and provide the following
 
 public:
 template<bool warn> struct Msg_;
 bool HasField(KeyCRef) const;
 std::string GetString(KeyCRef) const;
 void SetString(KeyCRef, std::string);
 typedef double ScalarType; // Or any suitable ScalarType

 protected:
 template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;
 template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
 */


template<typename Base> struct IO_ : Base {
    typedef Base Superclass;
//    Redeclare5Types(FromSuperclass,KeyType,KeyCRef,DiscreteType,ScalarType,ArrayOrdering);
    Redeclare5Types(Superclass,KeyType,KeyCRef,DiscreteType,ScalarType,ArrayOrdering);
    template<typename TC, size_t VD> using Array = typename Superclass::template Array<TC,VD>;
    template<size_t VD> using DimType = typename Superclass::template DimType<VD>;
    
    using Base::Base;
    
    typedef _Msg<false,IO_> Msg;
    typedef _Msg<true,IO_> WarnMsg;
	
    using Base::GetString;
    std::string GetString(KeyCRef, const std::string &, int=1) const;
    template<typename T> T Get(KeyCRef) const;
    template<typename T> T Get(KeyCRef, const T &, int=1) const;
    template<typename T> std::vector<T> GetVector(KeyCRef) const;
    template<typename T, size_t d> Array<T, d> GetArray(KeyCRef,
					ArrayOrdering=ArrayOrdering::Ignore) const;
    template<typename T> std::vector<DiscreteType> GetDimensions(KeyCRef,
					ArrayOrdering=ArrayOrdering::Ignore) const;
	int GetElementSize(KeyCRef, int) const;
    void SetHelp(KeyCRef, const std::string &);
    std::set<std::string> keyHelp;
    void UsageReport();

    // Output
    template<typename T> void Set(KeyCRef, const T &);
    template<typename T> void SetVector(KeyCRef, const std::vector<T> &);
    template<typename T, size_t d> void SetArray(KeyCRef, const Array<T, d> &,
					ArrayOrdering=ArrayOrdering::Ignore);
protected:
    template<typename V> static V TransposeDims(const V &);
    template<typename V> static V ReverseDims(const V &);
    template<typename T, size_t d, typename TI = const T *> struct TransposeVals;
    template<typename T, size_t d, typename TI = const T *> struct ReverseVals;
    template<typename T, size_t d, typename TI = const T *> struct TransposeReverseVals; // Transpose(Reverse())
    template<typename T, size_t d, typename TI = const T *> struct ReverseTransposeVals;
    
    template<typename T, size_t d> void Set(KeyCRef, DimType<d>, const T*);
    std::vector<std::string> unusedHelp;
};



#include "IO.hxx"

#endif /* IO_h */
