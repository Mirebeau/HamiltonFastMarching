// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef Curvature2Specializations_h
#define Curvature2Specializations_h

// This macro activates variants of the models which only allow a specific curvature sign
#ifndef convex_curvature_macro
#define convex_curvature_macro 0
#endif

#include "CommonTraits.h"

// --------------- R2S1 Traits ---------------

struct TraitsR2S1 : TraitsBase<3> {
    typedef EulerianDifference<OffsetType,ScalarType,1> DifferenceType;

    static const DiscreteType nStencilDependencies=1;
    constexpr static std::array<DiscreteType, nStencilDependencies>
    stencilDependencies = {{2}};
    
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Periodic}};
    typedef PeriodicGrid<TraitsR2S1> DomainType;
};
constexpr decltype(TraitsR2S1::stencilDependencies) TraitsR2S1::stencilDependencies;
constexpr decltype(TraitsR2S1::boundaryConditions) TraitsR2S1::boundaryConditions;


// --------------- 2D ReedsShepp model -------------
struct TraitsReedsShepp2 : TraitsR2S1 {
    typedef EulerianStencil<DifferenceType,4> StencilType;
};

struct StencilReedsShepp2 final
: HamiltonFastMarching<TraitsReedsShepp2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsShepp2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType xi=1; // xi is the typical curvature radius
    bool projective=false;
    
    typedef Traits::BasisReduction<2> ReductionType;
    Voronoi1Vec<ReductionType> reduc;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = index[2]*param.dependScale;
        const ReductionType::VectorType v{cos(theta),sin(theta)};
        
        auto & symmetric = stencil.symmetric[0];
        reduc(&symmetric[0],v/param.gridScale);
//        Voronoi1Mat<ReductionType>(&symmetric[0],
//        (Sym::RankOneTensor(v)*(1.-square(eps)) + Sym::Identity() * square(eps))/square(param.gridScale) );
        
        symmetric[3].offset = OffsetType{0,0,1};
        symmetric[3].baseWeight = 1./square(xi*param.dependScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc.eps = that->io.template Get<ScalarType>("eps",reduc.eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        projective=(bool)that->io.template Get<ScalarType>("projective",(ScalarType)projective);
        param.Setup(that,(projective ? mathPi : 2*mathPi)/dims.back());}
};

// --------------- 2D ReedsSheppForward model -------------
struct TraitsReedsSheppForward2 : TraitsR2S1 {
#if convex_curvature_macro
    typedef EulerianStencil<DifferenceType,0,4> StencilType;
#else
	typedef EulerianStencil<DifferenceType,1,3> StencilType;
#endif
};

struct StencilReedsSheppForward2 final
: HamiltonFastMarching<TraitsReedsSheppForward2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForward2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType xi=1; // xi is the typical curvature radius
    
    typedef Traits::BasisReduction<2> ReductionType;
    Voronoi1Vec<ReductionType> reduc;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = index[2]*param.dependScale;
        const typename ReductionType::VectorType v{cos(theta),sin(theta)};
        
        auto &forward = stencil.forward[0];
        reduc(&forward[0],v/param.gridScale);
        
#if convex_curvature_macro
		auto &angular = stencil.forward[0][3];
#else
		auto &angular = stencil.symmetric[0][0];
#endif
        angular.offset = OffsetType{0,0,-1};
        angular.baseWeight = 1./square(xi*param.dependScale);
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc.eps=that->io.template Get<ScalarType>("eps",reduc.eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/dims.back());}
};



// --------------- 2D Dubins model ---------------

struct TraitsDubins2 : TraitsR2S1 {
    typedef EulerianStencil<DifferenceType,0,6,2> StencilType;
};

struct StencilDubins2 final
: HamiltonFastMarching<TraitsDubins2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDubins2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType xi=1.;
    
    typedef Traits::BasisReduction<3> ReductionType;
    Voronoi1Vec<ReductionType> reduc;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = index[2]*param.dependScale;
        const VectorType
        v{cos(theta)/param.gridScale,sin(theta)/param.gridScale,1./(xi*param.dependScale)};
        reduc(&stencil.forward[0][0],v);
        
#if convex_curvature_macro
		const VectorType v0{cos(theta)/param.gridScale,sin(theta)/param.gridScale,0.};
        reduc(&stencil.forward[1][0],v0);
#else
        for(int i=0; i<6; ++i) {
            stencil.forward[1][i]=stencil.forward[0][i];
            stencil.forward[1][i].offset[2]*=-1;}
#endif
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc.eps=that->io.template Get<ScalarType>("eps",reduc.eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/dims.back());}
};
// ------------- 2D Elastica model ----------------

// Number of Fejer numerical integration points.

template<int nFejer>
struct TraitsElastica2 : TraitsR2S1 {
    static const DiscreteType nForward=6*nFejer - 3*((nFejer%2)==1);
    typedef EulerianStencil<DifferenceType,0,nForward> StencilType;
};

template<int nFejer>
struct StencilElastica2 final
: HamiltonFastMarching<TraitsElastica2<nFejer> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsElastica2<nFejer> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare7Types(Superclass,Traits,ScalarType,IndexType,VectorType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(Traits,mathPi)
    
    typename HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1.;
    static const std::array<ScalarType,nFejer> fejerWeights;
    
    template<size_t n> using BasisReduction = typename Traits::template BasisReduction<n>;
    Voronoi1Vec<BasisReduction<2> > reduc2;
    Voronoi1Vec<BasisReduction<3> > reduc3;

    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = (2.*mathPi*index[2])/this->dims[2];
        auto & forward = stencil.forward[0];
        for(int l=0; l<nFejer/2; ++l){
            const ScalarType phi = mathPi*(l+0.5)/nFejer;
            const VectorType v{
                sin(phi)*cos(theta)/this->param.gridScale,
                sin(phi)*sin(theta)/this->param.gridScale,
                cos(phi)/(xi*this->param.dependScale)
            };

            const int s0=2*l,s1=2*l+1;
            reduc3(&forward[6*s0], v);
            for(int i=0; i<6; ++i) {
                forward[6*s0+i].baseWeight*=fejerWeights[l];
                forward[6*s1+i] = forward[6*s0+i];
                forward[6*s1+i].offset[2]*=-1;
#if convex_curvature_macro
				forward[6*s1+i].baseWeight=0;
#endif
            } // for i
        } // for l
        if(nFejer%2==1){
            typedef typename Traits::template BasisReduction<2> ReductionType;
            typename ReductionType::VectorType
            v{cos(theta)/this->param.gridScale,sin(theta)/this->param.gridScale};
            const int s = nFejer-1;
            reduc2(&forward[6*s], v);
//            Voronoi1Vec<ReductionType>(&forward[6*s], v, eps);
            for(int i=0; i<3; ++i){
                // Multiplicative factor Scalar(sqrt(4/3.)) omitted 
                forward[6*s+i].baseWeight*=fejerWeights[nFejer/2];
#if convex_curvature_macro
				forward[6*s+i].baseWeight/=2;
#endif
				
			}
        } // if odd
		
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc3.eps=that->io.template Get<ScalarType>("eps",reduc3.eps);
        reduc2.eps=reduc3.eps;
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/this->dims.back());}
};

template<> const std::array<double, 1> StencilElastica2<1>::fejerWeights = {{2.}};
template<> const std::array<double, 2> StencilElastica2<2>::fejerWeights = {{1.,1.}};
template<> const std::array<double, 3> StencilElastica2<3>::fejerWeights =
{{0.444444, 1.11111, 0.444444}};
template<> const std::array<double, 4> StencilElastica2<4>::fejerWeights =
{{0.264298, 0.735702, 0.735702, 0.264298}};
template<> const std::array<double, 5> StencilElastica2<5>::fejerWeights =
{{0.167781, 0.525552, 0.613333, 0.525552, 0.167781}};
template<> const std::array<double, 6> StencilElastica2<6>::fejerWeights =
{{0.118661, 0.377778, 0.503561, 0.503561, 0.377778, 0.118661}};
template<> const std::array<double, 7> StencilElastica2<7>::fejerWeights =
{{0.0867162, 0.287831, 0.398242, 0.454422, 0.398242, 0.287831, 0.0867162}};
template<> const std::array<double, 8> StencilElastica2<8>::fejerWeights =
{{0.0669829, 0.222988, 0.324153, 0.385877, 0.385877, 0.324153, 0.222988, 0.0669829}};
template<> const std::array<double, 9> StencilElastica2<9>::fejerWeights =
{{0.0527366, 0.179189, 0.264037, 0.330845, 0.346384, 0.330845, 0.264037, 0.179189, 0.0527366}};

#endif /* Curvature2_h */
