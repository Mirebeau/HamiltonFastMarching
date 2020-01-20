// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_VectorType_h
#define AmongConvex2_VectorType_h

#include "PointType.h"
#include "FriendOperators.h"
#include "math.h"

namespace LinearAlgebra {
    
template <typename TComponent, size_t VDimension>
struct Vector :
public PointBase<TComponent, VDimension>,
affine_space<Point<TComponent, VDimension>, Vector<TComponent, VDimension> , TComponent>
{
    typedef TComponent ComponentType;
    static const size_t Dimension = VDimension;
    
    typedef Point<ComponentType, Dimension> PointType;
    typedef PointBase<ComponentType, Dimension> PointBaseType;
    
    typedef typename PointBaseType::DataType DataType;
    template<typename ...T,typename dummy = typename std::enable_if<sizeof...(T)==Dimension>::type >
    constexpr Vector(T... t):PointBaseType{t...}{};

    Vector(){};
    
    // Linear algebra
    Vector & operator+=(const Vector & u){
		for(int i=0; i<Dimension; ++i) this->operator[](i)+=u[i]; return *this;}
    Vector & operator-=(const Vector & u){
		for(int i=0; i<Dimension; ++i) this->operator[](i)-=u[i]; return *this;}
    Vector & operator*=(const ComponentType & a){
		for(int i=0; i<Dimension; ++i) this->operator[](i)*=a; return *this;}
    Vector & operator/=(const ComponentType & a){
		if constexpr(std::is_same_v<decltype(a!=ComponentType(0.)),bool>) assert(a != ComponentType(0.));
        const ComponentType b=ComponentType(1)/a;
        for(int i=0; i<Dimension; ++i) this->operator[](i)*=b;
        return *this;
    }
    
    Vector operator -() const {Vector u; for(int i=0; i<Dimension; ++i) u[i]=-this->operator[](i); return u;}
    
    friend PointType & operator+=(PointType & p, const Vector & u){for(int i=0; i<Dimension; ++i) p[i]+=u[i]; return p;}
    friend PointType & operator-=(PointType & p, const Vector & u){for(int i=0; i<Dimension; ++i) p[i]-=u[i]; return p;}
    
    // Geometry
	template<typename T, typename S=algebra_t<ComponentType, T> >
    S ScalarProduct(const Vector<T,Dimension> & u) const {S s(0.);
		for(size_t i=0; i<Dimension; ++i) s += this->operator[](i)*u[i]; return s;}
    ComponentType SquaredNorm() const {return ScalarProduct(*this);}
	ComponentType Norm() const { assert(!std::is_integral<ComponentType>::value); using std::sqrt;  return sqrt(SquaredNorm()); }
    bool IsNull() const {for(int i=0;i<Dimension;++i) if((*this)[i]!=ComponentType(0.)) return false; return true;}
    Vector Normalized() const;
	
	Vector ComponentWiseProduct(const Vector & other) const {
		Vector res; for(int i=0; i<Dimension; ++i) res[i]=this->operator[](i)*other[i]; return res;}

    // Constructors
    static Vector FromOrigin(const PointBaseType & p){return Vector(p);}
    static Vector Constant(ComponentType c){return PointBaseType::Constant(c);}
    static Vector RandomUnit();
    template<typename Component2>
    static Vector CastCoordinates(const Vector<Component2,Dimension> & u){return Vector(PointBaseType::CastCoordinates(u));}
	static constexpr Vector CanonicalBasis(int i){assert(0<=i && i<Dimension);
		Vector v; v.fill(0); v[i]=1; return v;}
protected:
    Vector(const PointBaseType & p):PointBaseType(p){};
};

    
template<typename TC, size_t VD>
Vector<TC,VD>  operator- (const Point<TC,VD> & p, const Point<TC,VD> & q)
{
    Vector<TC,VD> u;
    for(int i=0; i<VD; ++i)
        u[i]=p[i]-q[i];
    return u;
}

    
// Determinant, in dimension 2,3.
template<typename T1,typename T2,typename S = algebra_t<T1, T2> >
S Determinant(const Vector<T1,2> & u, const Vector<T2,2> & v){return u[0]*v[1]-u[1]*v[0];}

template<typename TComponent>
TComponent
Determinant(const Vector<TComponent,3> & u, const Vector<TComponent,3> & v, const Vector<TComponent,3> & w){
    TComponent ans=0;
    for(int i=0; i<3; ++i)
        ans+= u[i]*v[(i+1)%3]*w[(i+2)%3] - u[(i+2)%3]*v[(i+1)%3]*w[i];
    return ans;
}
    
// Perpendicular and cross product
template<typename TComponent>
Vector<TComponent, 2>
Perp(const Vector<TComponent, 2> & v){
    return {-v[1],v[0]};
}
    
template<typename TComponent>
Vector<TComponent, 3>
Cross(const Vector<TComponent, 3> & u, const Vector<TComponent, 3> & v){
    Vector<TComponent, 3> result;
    for(int i=0; i<3; ++i) {
        int j=(i+1)%3, k=(i+2)%3;
        result[i] = u[j]*v[k]-u[k]*v[j];
    }
    return result;
}
    

    
template<typename TComponent>
double
AngleWithRespectToNegativeAxis(const Vector<TComponent, 2> & u){return atan2(u[1],u[0]);}

#include "Implementation/VectorType.hxx"

} // end of namespace LinearAlgebra

#endif
