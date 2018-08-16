// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Curvature3_h
#define Curvature3_h

#include "CommonTraits.h"

// ------------- R3S2 Traits ---------------

struct TraitsR3S2 : TraitsBase<5> {
    typedef Difference<1> DifferenceType;
    static const DiscreteType nStencilDependencies=2;
    constexpr static std::array<DiscreteType, nStencilDependencies>
    stencilDependencies = {{3,4}};
    
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed,
        Boundary::Sphere2_0, Boundary::Sphere2_1}};
    
    typedef PeriodicGrid<TraitsR3S2> DomainType;
};
constexpr decltype(TraitsR3S2::stencilDependencies) TraitsR3S2::stencilDependencies;
constexpr decltype(TraitsR3S2::boundaryConditions) TraitsR3S2::boundaryConditions;

// ---------- 3D ReedsSheppForward. Includes the dual torsion-like model. ----------

struct TraitsReedsShepp3 : TraitsR3S2 {
    static const DiscreteType nSymmetric=8;
};

struct StencilReedsShepp3
: HamiltonFastMarching<TraitsReedsShepp3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsShepp3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1;
    bool projective=false, dual=false;
    
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType
        theta=(index[3]+0.5)*param.dependScale,
        phi=index[4]*param.dependScale;
        typedef Traits::BasisReduction<3> ReductionType;
        typedef typename ReductionType::VectorType VectorType;
        typedef typename ReductionType::SymmetricMatrixType Sym;
        const VectorType v =
        (VectorType{cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi)});
        const Sym m = dual ?
        (square(eps)-1.)*Sym::RankOneTensor(v)+ Sym::Identity() :
        (1.-square(eps))*Sym::RankOneTensor(v)+square(eps)*Sym::Identity();
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m/square(param.gridScale));
        //        Voronoi1Vec<ReductionType>(&stencil.symmetric[0], v/param.gridScale, eps);

        stencil.symmetric[6].offset = OffsetType{0,0,0,1,0};
        stencil.symmetric[6].baseWeight = 1./square(xi*param.dependScale);
        
        stencil.symmetric[7].offset = OffsetType{0,0,0,0,1};
        stencil.symmetric[7].baseWeight = 1./square(sin(theta)*xi*param.dependScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {auto & io = that->io;
        Superclass::Setup(that);
        eps=io.template Get<ScalarType>("eps",eps);
        xi=io.template Get<ScalarType>("xi",xi);
        dual = (bool) io.template Get<ScalarType>("dual",dual);
        projective= (bool) io.template Get<ScalarType>("projective",projective);
        param.Setup(that,2*mathPi/dims.back()); 
        
        if(projective && dims[4]!=4*dims[3])
            ExceptionMacro("ReedsShepp3 setup error : inconsistent bundle dimensions (projective : nPhi=4*nTheta).\n");
        if(!projective && dims[4]!=2*dims[3])
            ExceptionMacro("ReedsShepp3 setup error : inconsistent bundle dimensions (expected : nPhi=2*nTheta).\n");
    }
};

// ---------- 3D ReedsSheppForward model ----------

struct TraitsReedsSheppForward3 : TraitsR3S2 {
    static const DiscreteType nForward=6, nSymmetric=2;
};

struct StencilReedsSheppForward3
: HamiltonFastMarching<TraitsReedsSheppForward3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForward3> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType
        theta = mathPi*(index[3]+0.5)/dims[3],
        phi = (2.*mathPi*index[4])/dims[4];
        typedef Traits::BasisReduction<3> ReductionType;
        typedef typename ReductionType::VectorType VectorType;
        const VectorType v =
        (VectorType{cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi)})/param.gridScale;
        Voronoi1Vec<ReductionType>(&stencil.forward[0], v, eps);
        
        stencil.symmetric[0].offset = OffsetType{0,0,0,1,0};
        stencil.symmetric[0].baseWeight = 1./square(xi*param.dependScale);
        
        stencil.symmetric[1].offset = OffsetType{0,0,0,0,1};
        stencil.symmetric[1].baseWeight = 1./square(sin(theta)*xi*param.dependScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {auto & io = that->io;
        Superclass::Setup(that);
        eps=io.template Get<ScalarType>("eps",eps);
        xi=io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/dims.back());
        if(dims[4]!=2*dims[3])
            ExceptionMacro("ReedsSheppForward3 setup error : inconsistent bundle dimensions (expected : nPhi=2*nTheta).\n");
    }
};

#endif /* Curvature3_h */
