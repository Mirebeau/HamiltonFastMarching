// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef JMMCppLibs_ExportMacros_h
#define JMMCppLibs_ExportMacros_h

/*

Text based export of variables, arrays, images. Mathematica friendly output.
Also defines an << operator for printing pairs in the form {first,second}.
 
std::cout << "{"
    ExportVarArrow(var)
    ExportNamedVarArrow(var,fancyName)
    ExportArrayArrow(array)
    ExportArrayRecursiveArrow(arrayOfArrays,depth)
    ExportNamedVarArrow(MakeImageRegion(begin,end,dimension,image),fancyName)
    << "}";
 
 ExportArrayRecursiveArrow is intended for arrays of arrays of ..., up to specified depth.
 
 MakeImageRegion is intended for multidimensional arrays, accessed via multi-indices.
 it = (it1,it2,...,itd); value = image(it);
 
 */

#include <iostream>
#include <iterator>

// ----- In class << operator declaration ------

#define PrintSelfMacro(type) \
void PrintSelf(std::ostream &) const; \
friend std::ostream & operator << (std::ostream & os, const type & a){a.PrintSelf(os); return os;}

// --------- << operator for printing pairs -------

template<typename S, typename T>
std::ostream & operator << (std::ostream & os, const std::pair<S,T> & p){
    return os << "{" << p.first << "," << p.second << "}";}

// --------- Export single variable ---------

#define ExportVarArrow(name) \
<< '"' << #name << '"' << " -> "  <<  name  << ","

#define ExportEnumArrow(name) \
<< '"' << #name << '"' << " -> "  <<  enumToString(name)  << ","

#define ExportNamedVarArrow(name,var) \
<< '"' << name << '"' << " -> "  <<  var  << ","

// --------- Export array of values -----------

template<int depth, typename DataType>
struct _ExportArrayContainer {
	const DataType & data;};

template<int depth, typename DataType>
_ExportArrayContainer<depth,DataType> _Make_ExportArrayContainer(const DataType & data){
    return _ExportArrayContainer<depth,DataType>{data};}

template<int depth, typename DataType>
std::ostream & operator << (std::ostream & os, _ExportArrayContainer<depth,DataType> container){
    auto it=std::begin(container.data);
	os << "{";
    if(it!=std::end(container.data))
		while(true){
            os << _Make_ExportArrayContainer<depth-1>(*it);
			++it;
			if(it!=std::end(container.data)) os << ",";
			else break;
		}
	os << "}";
	return os;
}

template<typename DataType>
std::ostream & operator << (std::ostream & os, _ExportArrayContainer<0,DataType> container){
    auto it=std::begin(container.data);
    os << "{";
    if(it!=std::end(container.data))
        while(true){
            os << *it;
            ++it;
            if(it!=std::end(container.data)) os << ",";
            else break;
        }
    os << "}";
    return os;
}

template<typename TFirst,typename TSecond>
std::ostream & operator << (std::ostream & os, _ExportArrayContainer<0, std::pair<TFirst,TSecond> > container){
    return os << "{" << container.data.first << "," << container.data.second << "}";
}

#define ExportArrayArrow(name) \
<< '"' << #name << '"' << " -> " << _Make_ExportArrayContainer<0>(name) << ","

#define ExportArrayRecursiveArrow(name,depth) \
<< '"' << #name << '"' << " -> " << _Make_ExportArrayContainer<depth>(name) << ","
 
// ------- Export image -------

template<typename ImageType, typename RegionIteratorType>
struct _ImageRegion {
	RegionIteratorType begin, end;
	int Dimension;
	ImageType im;
};

template<typename ImageType, typename RegionIteratorType>
_ImageRegion<ImageType,RegionIteratorType> MakeImageRegion(
		RegionIteratorType begin, RegionIteratorType end, int Dimension, ImageType im){
	return _ImageRegion<ImageType,RegionIteratorType>{begin,end,Dimension,im};}

template<typename ImageType, typename RegionIteratorType>
std::ostream & operator << (std::ostream & os, _ImageRegion<ImageType,RegionIteratorType> region){
	for(int i=0; i<region.Dimension; ++i) os << "{";
	for(int i=0; i<region.Dimension; ++i)  // empty array case
		if(region.begin[i]>=region.end[i]){
			for(int i=0; i<region.Dimension; ++i) os << "}";
			return os;}
	RegionIteratorType it = region.begin;
	for(int i=0; i<region.Dimension; ){
		os << region.im(it);
		for(i=0; i<region.Dimension; ++i){
			if(i>0) os << "}";
			if(it[i]+1<region.end[i]){
				++it[i];
				if(it[i]<region.end[i])
					os << ",";
				for(int j=0; j<i; ++j){
					it[j]=region.begin[j];
					os << "{";
				}
				break;
			}
		}
	}
	return os << "}";
}

// --------- Print a range of values -------

template <typename ForwardIterator>
void print_range(std::ostream & f, ForwardIterator first, ForwardIterator last) {
    f<<"{";
    if(first!=last) f<<*first++;
    while(first!=last) f << "," << *first++;
    f<<"}";
}


#endif
