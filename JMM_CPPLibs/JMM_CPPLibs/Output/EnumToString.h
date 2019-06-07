// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

// This is a modified version of the code appearing in the thread
// http://codereview.stackexchange.com/questions/14309/conversion-between-enum-and-string-in-c-class-header

#ifndef EnumToString_h
#define EnumToString_h

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <sstream>

/* //Usage
 struct Y {
 enum class X {Hi,Lo};
 };
 template<> char const* enumStrings<Y::X>::data[] = {"Hi", "Lo"};
 */

// This is the type that will hold all the strings.
// Each enumerate type will declare its own specialization.
// Any enum that does not have a specialization will generate a compiler error
// indicating that there is no definition of this variable (as there should be
// be no definition of a generic version).
template<typename T>
struct enumStrings
{
    static char const* data[];
};

// This is a utility type.
// Creted automatically. Should not be used directly.
template<typename T>
struct enumRefHolder
{
    T& enumVal;
    enumRefHolder(T& enumVal): enumVal(enumVal) {}
};
template<typename T>
struct enumConstRefHolder
{
    T const& enumVal;
    enumConstRefHolder(T const& enumVal): enumVal(enumVal) {}
};

// The next too functions do the actual work of reading/writtin an
// enum as a string.
template<typename T>
std::ostream& operator<<(std::ostream& str, enumConstRefHolder<T> const& data)
{
    return str << enumStrings<T>::data[static_cast<int>(data.enumVal)];
}

template<typename T>
std::istream& operator>>(std::istream& str, enumRefHolder<T> const& data)
{
    std::string value;
    str >> value;
    
    static auto begin  = std::begin(enumStrings<T>::data);
    static auto end    = std::end(enumStrings<T>::data);
    
    auto find   = std::find(begin, end, value);
    if (find != end)    data.enumVal = static_cast<T>(std::distance(begin, find));
    else                data.enumVal = static_cast<T>(-1);
    return str;
}


// This is the public interface:
// use the ability of function to deduce their template type without
// being explicitly told to create the correct type of enumRefHolder<T>
template<typename T>
enumConstRefHolder<T>  enumToString(T const& e) {return enumConstRefHolder<T>(e);}

template<typename T>
enumRefHolder<T>       enumFromString(T& e)     {return enumRefHolder<T>(e);}

template<typename T>
std::string enumToRealString(T const& e){
    std::ostringstream oss;
    oss << enumToString(e);
    return oss.str();
}

template<typename T>
T enumFromString(std::string value){
    static auto begin  = std::begin(enumStrings<T>::data), end = std::end(enumStrings<T>::data);
    auto find   = std::find(begin, end, value);
    if (find != end) return static_cast<T>(std::distance(begin, find));
    else             return static_cast<T>(-1);    
}
#endif
