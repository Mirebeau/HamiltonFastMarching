// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_PointType_h
#define AmongConvex2_PointType_h

#include "PointBaseType.h"

namespace LinearAlgebra {
    
template<typename TComponent, size_t VDimension>
struct Point :
PointBase<TComponent,VDimension>
{
    typedef TComponent ComponentType;
    static const int Dimension = VDimension;
    typedef PointBase<TComponent,VDimension> PointBaseType;
    
    template<typename ...T>//,typename dummy = typename std::enable_if<sizeof...(T)==Dimension>::type >
    constexpr Point(T... t):PointBaseType{t...}{};

    Point(){};
    
    static Point FromOrigin(const PointBaseType & p){return Point(p);}
    static Point Constant(ComponentType c){return PointBaseType::Constant(c);}
    template<typename ComponentType2>
    static Point CastCoordinates(const Point<ComponentType2,Dimension> & p){
        return Point(PointBaseType::CastCoordinates(p));}
    
    template<size_t n>
    static Point
    Barycenter(const std::array<Point,n> & points, const std::array<TComponent,n> & weights) {
        assert(fabs(std::accumulate(weights.begin(), weights.end(), ComponentType(0)) - 1) < 1e-6);
        Point p;
        for(int i=0; i<Dimension; ++i){
            p[i]=0;
            for(int j=0; j<n; ++j)
                p[i]+=weights[j]*points[j][i];
        }
        return p;
    }
    
    static Point Barycenter(const Point & p, const Point & q, ComponentType t) {
        assert(0<=t && t<=1);
        Point b;
        for(int i=0; i<Dimension; ++i)
            b[i]=(1-t)*p[i]+t*q[i];
        return b;
    }
protected:
    Point(const PointBaseType & p): PointBaseType(p){};
};
    
}
#endif
