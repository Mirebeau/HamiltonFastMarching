// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef QuaternionSpecializations_h
#define QuaternionSpecializations_h

#include "Specializations/CommonTraits.h"

/* This file was intended for curvature penalized shortest path on S2.
 It does not work so weel because the combination of the singularity of the parametrization,
 which leads to unreasonnable condition numbers and stencil widths. 
 */


// ---------- SO3 Traits ----------
struct TraitsSO3 : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Sphere3_0, Boundary::Proj3_1, Boundary::Sphere3_2}};
};
constexpr decltype(TraitsSO3::boundaryConditions) TraitsSO3::boundaryConditions;


// Parametrization of SO3 by three angles.
// Regarded as a {point in S2, subriemannian tangent, sideways direction}
// Domain : t1 in [0,Pi/2], t2 in [0,Pi] and t3 in [0,2Pi].
/*
template<typename VectorType,
typename ScalarType=typename VectorType::ScalarType>
std::array<VectorType,3> SphericalPoint(ScalarType t1, ScalarType t2, ScalarType t3){
    const ScalarType
    c1=cos(t1), s1=sin(t1),
    c2=cos(t2), s2=sin(t2),
    c3=cos(t3), s3=sin(t3);
    const ScalarType
    a=c1, b=s1*c2, c=s1*s2*c3, d=s1*s2*s3;
    
    return {VectorType{square(a)+square(b)-square(c)-square(d),2*b*c-2*a*d,2*a*c+2*b*d},
        VectorType{2*b*c+2*a*d,square(a)-square(b)+square(c)-square(d),-2*a*b+2*c*d},
        VectorType{-2*a*c+2*b*d,2*a*b+2*c*d,square(a)-square(b)-square(c)+square(d)}};
}*/

// Vector fields, in the angular coordinates, for going {forward, rotating, stepping aside}.
template<typename VectorType,
typename ScalarType=typename VectorType::ScalarType>
std::array<VectorType,3> SphericalFields(ScalarType t1, ScalarType t2, ScalarType t3){
    const ScalarType
    c1=cos(t1), s1=sin(t1),
    c2=cos(t2), s2=sin(t2),
    c3=cos(t3), s3=sin(t3);
    
    const ScalarType
    m11 = -s1*c2, m12 = -c1*s2,
    m21 = -s1*s2, m22 = c1*c2,
    m31 = c1*c3,        m33 = -s1*s3,
    m41 = c1*s3,        m43 = s1*c3;
    
    const ScalarType d1=1, d2=square(c1), d3=square(s1);
    
    auto chgCoord = [m11,m12,m21,m22,m31,m33,m41,m43,d1,d2,d3]
    (const std::array<ScalarType,4> & v)->VectorType {
        return VectorType{(m11*v[0]+m21*v[1]+m31*v[2]+m41*v[3])/d1,
            (m12*v[0]+m22*v[1])/d2,
            (m33*v[2]+m43*v[3])/d3};};
    
    const ScalarType a=c1*c2, b=c1*s2, c=s1*c3, d=s1*s3;
    
/*  //Old, too much singular, parametrization
 const ScalarType
    m11=-s1,    m12=c1*c2,  m13=c1*s2*c3,   m14=c1*s2*s3,
                m22=-s1*s2, m23=s1*c2*c3,   m24=s1*c2*s3,
                            m33=-s1*s2*s3,  m34=s1*s2*c3;
    
    const ScalarType d1=1,d2=square(s1),d3=square(s1*s2);
    
    auto chgCoord = [m11,m12,m13,m14,m22,m23,m24,m33,m34,d1,d2,d3](std::array<ScalarType,4> v)->VectorType{
        return VectorType{
            (m11*v[0]+m12*v[1]+m13*v[2]+m14*v[3])/d1,
                     (m22*v[1]+m23*v[2]+m24*v[3])/d2,
                              (m33*v[2]+m34*v[3])/d3};
    };
    
    const ScalarType
    a=c1, b=s1*c2, c=s1*s2*c3, d=s1*s2*s3;*/
    
    return {{
        chgCoord({{d, c, -b, -a}}),
        chgCoord({{b, -a, d, -c}}),
        chgCoord({{-c, d, a, -b}}) }};
}

// ---------- ReedsShepp model on SO3 ---------

struct TraitsReedsSheppSO3 : TraitsSO3 {
    static const DiscreteType nSymmetric = 6;
};

struct StencilReedsSheppSO3
: HamiltonFastMarching<TraitsReedsSheppSO3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppSO3> HFM;
    typedef HamiltonFastMarching<TraitsReedsSheppSO3>::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1;
    typedef ScalarType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        const ScalarType h=param.gridScale;
        const auto v=SphericalFields<VectorType>((index[0]+0.5)*h,index[1]*h, index[2]*h);
        typedef Traits::BasisReduction<3> ReductionType;
        typedef ReductionType::SymmetricMatrixType Sym;
        const Sym diff = Sym::RankOneTensor(v[0])+Sym::RankOneTensor(v[1]/xi)+Sym::RankOneTensor(eps*v[2]);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], diff/square(h));
/*        const ScalarType theta=(index[0]+0.5)*h, phi=(index[1]+0.5)*h;
        assert(theta<mathPi/2 && phi<mathPi);
        stencil.symmetric[0].offset={1,0,0};
        stencil.symmetric[0].baseWeight=1;
        stencil.symmetric[1].offset={0,1,0};
        stencil.symmetric[1].baseWeight=1./square(sin(theta));
        stencil.symmetric[2].offset={0,0,1};
        stencil.symmetric[2].baseWeight=1./square(sin(theta)*sin(phi));
        for(int i=3; i<6; ++i){
            stencil.symmetric[i].offset.fill(0);
            stencil.symmetric[i].baseWeight=0;
        }*/
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        pDualMetric = that->template GetField<ScalarType>("dualMetric");
        const int n=dims[0];
        if(dims[1]!=2*n || dims[2]!=4*n)
            ExceptionMacro("Inconsistent domain dims: should be (n/2,n,2n)");
        const ScalarType h=Traits::mathPi/(2.*n);
        param.gridScale = h;
        param.origin[1]=-h/2; param.origin[2]=-h/2;
    }
};


// ---------------  ReedsSheppForward model on SO3 -------------
struct TraitsReedsSheppForwardSO3 : TraitsSO3 {
    static const DiscreteType nForward=6, nSymmetric=6;
};

struct StencilReedsSheppForwardSO3
: HamiltonFastMarching<TraitsReedsSheppForwardSO3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForwardSO3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1; // xi is the typical curvature radius
    typedef ScalarType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType h=param.gridScale;
        const auto & v = SphericalFields<VectorType>((index[0]+0.5)*h,(index[1]+0.5)*h, (index[2]+0.5)*h); // TODO : what origin ?
        
        typedef Traits::BasisReduction<3> ReductionType;
        Voronoi1Vec<ReductionType>(&stencil.forward[0], v[0]/h, eps);
        Voronoi1Vec<ReductionType>(&stencil.symmetric[0], v[1]/(xi*h), eps);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        pDualMetric = that->template GetField<ScalarType>("dualMetric");
        const int n=dims[0];
        if(dims[1]!=2*n || dims[2]!=4*n)
            ExceptionMacro("Inconsistent domain dims: should be (n/2,n,2n)");
        const ScalarType h=Traits::mathPi/(2.*n);
        param.gridScale = h;
        param.origin[1]=-h/2; param.origin[2]=-h/2;
    }
};


// ---------------  Dubins model on SO3 -------------

/*
struct TraitsDubinsSO3 : TraitsSO3 {
    static const DiscreteType nMax=6, nMaxForward=6;
};

struct StencilDubinsSO3
: HamiltonFastMarching<TraitsReedsSheppForwardSO3>::StencilDataType {
    ParamType param;
    ScalarType eps=0.1, xi=1; // xi is the typical curvature radius
    typedef ScalarType MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType h=param.gridScale;
        const auto & v = SphericalFields<VectorType>((index[0]+0.5)*h,index[1]*h, index[2]*h);
        
        const auto
        w1 = (v[0]+v[1]/xi)/param.gridScale,
        w2 = (v[0]-v[1]/xi)/param.gridScale;
        
        typedef Traits::BasisReduction<3> ReductionType;
        Voronoi1Vec<ReductionType>(&stencil.forward[0], w1, eps);
        Voronoi1Vec<ReductionType>(&stencil.symmetric[0], w2, eps);
    }
};*/

// ------ TODO : Euler elastica ? Rolling ball ? --------


#endif /* QuaternionSpecializations_h */
