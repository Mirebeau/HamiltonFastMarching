//
//  RollingBall.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 15/02/2018.
//

#ifndef RollingBall_h
#define RollingBall_h

/*
 In this notebook we make a new attempt at problems solving problems involving the which configuration space of the ball S^2.
 
 A previous experiment involved the group SO3, and its parametrization by quaternions, aka S^4. This turned out to be a failure, due to the difficulty to parametrize S^4 with a cartesian grid.
 
 In this new attempt, we describe the sphere S2 using the usual coordinate system (cos theta, sin theta cos phi, sin theta sin phi). A tangent circle of orientations is put on top of each point, with the orientation origin aligned with the meridian defined by phi. (Which is not a geodesic.)

 At each point, we thus have three vector fields of interest, associated with the angle theta, phi, and psi. A suitable boundary condition is needed for psi
 
 */

// TODO : Allow decoupling the stencil structure from the speed input structure.

#include "Base/HamiltonFastMarching.h"
#include "Specializations/CommonTraits.h"

// ---------- Reeds-Shepp (projective) on the sphere -------

// Begin with a projective orientation space, which makes sense for the Reeds shepp model. Hence no need for weird boundary conditions
struct TraitsReedsSheppS2 : TraitsBase<3> {
    typedef Difference<1> DifferenceType;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Sphere2_0, Boundary::Sphere2_1, Boundary::Periodic}};
    
    static const DiscreteType nStencilDependencies=2;
    constexpr static std::array<DiscreteType, nStencilDependencies>
    stencilDependencies = {{0,2}};
    // TODO : allow 'stencil dependencies' to be distinct from 'speed dependencies' (?)
    
    static const DiscreteType nSymmetric = 3+1;
};
// Linker wants the following line for some obscure reason.
constexpr const decltype(TraitsReedsSheppS2::boundaryConditions) TraitsReedsSheppS2::boundaryConditions;

struct StencilReedsSheppS2
: HamiltonFastMarching<TraitsReedsSheppS2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppS2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::_ParamDefault<2> param; // Diagonal parametrization
    ScalarType eps=0.1, xi=1; // xi is the typical curvature radius
    bool projective=false;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const auto & scales = param.gridScales;
        const ScalarType theta = (index[0]+0.5)*scales[0];
//        const ScalarType phi = index[1]*scales[1];
        const ScalarType psi = index[2]*scales[2];

        typedef Traits::BasisReduction<3> ReductionType;
        const ReductionType::VectorType // Issue : is this primal or dual ??
        vTheta = {1,0,0},
        vPhi = {0,1/sin(theta),-cos(theta)/sin(theta)},
        vPsi = {0,0,xi};
        
        const ReductionType::VectorType
        u = cos(psi)*vTheta+sin(psi)*vPhi,
        v = -sin(psi)*vTheta+cos(psi)*vPhi;
        
        typedef ReductionType::SymmetricMatrixType Sym;
        // Define the metric
        Sym m = Sym::RankOneTensor(u) + Sym::RankOneTensor(v)/square(eps) + Sym::RankOneTensor(vPsi);
        for(int i=0; i<=Dimension; ++i){
            for(int j=0; j<=i; ++j){
                m(i,j)*=scales[i]*scales[j];}}
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],m.Inverse());
        
        /*
        typedef Traits::BasisReduction<2> ReductionType;
        const ReductionType::VectorType v{cos(psi)/scales[0],sin(psi)/(sin(theta)*scales[1])};
        typedef ReductionType::SymmetricMatrixType Sym;
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],
            (Sym::RankOneTensor(v)*(1.-square(eps)) + Sym::Identity() * square(eps)) );
        
        // Motion in the angular direction
        stencil.symmetric[3].offset = OffsetType{0,0,1};
        stencil.symmetric[3].baseWeight = 1./square(xi*param.gridScales[2]);*/
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.gridScales = {mathPi/dims[0], 2.*mathPi/dims[1], mathPi/dims[2]};
        param.origin = {0.,-param.gridScales[1]/2,-param.gridScales[2]/2};
    }
};

// ---------------- Dubins on the sphere -------------
struct TraitsDubinsS2 : TraitsBase<3> {
    typedef Difference<1> DifferenceType;
    constexpr static const std::array<Boundary, Dimension>  boundaryConditions =
    {{Boundary::Sphere2_0, Boundary::Sphere2_1, Boundary::Sphere2_Hopf}};

    static const DiscreteType nStencilDependencies=2;
    constexpr static std::array<DiscreteType, nStencilDependencies>
    stencilDependencies = {{0,2}};
    
    static const DiscreteType nMax = 2, nMaxForward = 6;
};
constexpr const decltype(TraitsDubinsS2::boundaryConditions) TraitsDubinsS2::boundaryConditions;

struct StencilDubinsS2
: HamiltonFastMarching<TraitsDubinsS2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDubinsS2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::_ParamDefault<2> param; // Diagonal parametrization
    ScalarType eps=0.1, xi=1; // xi is the typical curvature radius
    bool projective=false;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const auto & scales = param.gridScales;
        const ScalarType theta = (index[0]+0.5)*scales[0];
        //        const ScalarType phi = index[1]*scales[1];
        const ScalarType psi = index[2]*scales[2];
        
        typedef Traits::BasisReduction<3> ReductionType;
        for(int i=0; i<2; ++i){
            ReductionType::VectorType
            v{cos(psi),sin(psi)/sin(theta),-sin(psi)*cos(theta)/sin(theta)+(i==0 ? -1 : 1)/xi};
            for(int j=0; j<Dimension; ++j) v[j]/=scales[j];
            Voronoi1Vec<ReductionType>(&stencil.maxForward[i][0],v,eps);
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.gridScales = {mathPi/dims[0], 2.*mathPi/dims[1], 2.*mathPi/dims[2]};
        param.origin = {0.,-param.gridScales[1]/2,-param.gridScales[2]/2};
    }
};

// ----------------- Rolling ball on the plane -------
struct TraitsRollingBall : TraitsBase<5> {
    typedef  Difference<0> DifferenceType;
    constexpr static const std::array<Boundary,Dimension> boundaryConditions = {{
        Boundary::Closed, Boundary::Closed, // x,y
        Boundary::Sphere2_0, Boundary::Sphere2_1, Boundary::Sphere2_Hopf //theta,phi,psi
    }};
    
    static const DiscreteType nStencilDependencies = 3;
    constexpr static std::array<DiscreteType,nStencilDependencies>
    stencilDependencies = {{2,3,4}};
    
    static const DiscreteType nSymmetric = 15;
};
constexpr const decltype(TraitsRollingBall::boundaryConditions) TraitsRollingBall::boundaryConditions;

struct StencilRolling2
: HamiltonFastMarching<TraitsRollingBall>::StencilDataType {
    typedef HamiltonFastMarching<TraitsRollingBall> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::_ParamDefault<2> param; // Diagonal parametrization
    ScalarType eps=0.1, xi=1;
    bool projective=false;
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const auto & scales = param.gridScales;
        const ScalarType theta = (index[2]+0.5)*scales[2];
        //        const ScalarType phi = index[3]*scales[3];
        const ScalarType psi = index[4]*scales[4];
        
        typedef Traits::BasisReduction<5> ReductionType;
        const ReductionType::VectorType
        // roll in the direction theta
        v0{cos(psi),sin(psi),1,0,0},
        // roll in the direction phi
        v1{-sin(psi),cos(psi),0,1/sin(theta),-cos(theta)/sin(theta)},
        // roll in the direction psi (vertical axis)
        v2{0,0,0,0,xi};
        

        
        // Motion in the angular direction
        
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        eps=that->io.template Get<ScalarType>("eps",eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
//        param.gridScales = {mathPi/dims[0], 2.*mathPi/dims[1], 2.*mathPi/dims[2]};
//        param.origin = {0.,-param.gridScales[1]/2,-param.gridScales[2]/2};
    }
};

#endif /* RollingBall_h */
