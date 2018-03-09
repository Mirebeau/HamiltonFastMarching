// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef AsymmetricQuadratic_h
#define AsymmetricQuadratic_h

#include "Specializations/CommonTraits.h"
#include "LinearAlgebra/VectorPairType.h"

/*******
 This file implements hamiltonians of 'asymmetric quadratic' type. That is of the form
 2 H(g) := <g,D g> + <eta, g>_+^2
 
 The standard quadratic term <g,D g> is discretized in the usual Riemannian manner.
 The positive linear part <eta, g>_+^2 is discretized using a relaxation parameter eps.
 
 Note : in two dimensions, a consistent and parameterless discretization
 is proposed in the ITKFM, based on an entirely distinct approach.
 (A semi-Lagrangian discretization using the FM-ASR algorithm.)
 
 In practice, adequately defining the hamiltonian may be non-intuitive.
 We thus provide two helper methods.
 - The option to dualize, that is to define the Legendre-Fenchel dual
 2 L(v) := <v, M v> + <omega, v>_+^2.
 Algebraic formulas relate (M,eta) with (D,v).
 - The option to use half-disk models, suitably relaxed. That is
 L(v) := a^2 <omega,v>_+^2 + infty <omega,v>_-^2 + b^2 |P_eta(v)|^2.
 where eta is a unit vector. We denoted by P_eta the orthogonal projector onto eta^perp.
 A relaxation is also needed for this term.
 
 Input format:
 - Standard and Dual case:
 % Symmetric tensor components (D)
 xx,
 xy, yy,
 xz, yz, zz,
 % Followed by the symmetric vector components (eta)
 x, y, z
 
 % A single relaxation parameter eps.
 
 - Half-Disk case
 % Vector components of eta/a (foward speed)
 x, y, z
 % Scalar field s = a/b (sideways speed ratio)
 s

 Relaxation parameters are:
 - a scalar eps, for how badly is backwards motion penalized. (relaxation of |P_eta(v)|^2
 - a scalar epsForward, for how accurately is forward speed implemented.
 As usual, small parameters yield large stencils.
 
********/

namespace AsymmetricNorm {
// TODO a distinct grid scale for the radius dimension in AsymmetricQuadratic3p1.

// Utility functions, to convert between the different input types

// In place dualization of (M,v) into (D,w)
template<typename SymmetricMatrixType, typename VectorType>
void Dualize(SymmetricMatrixType & M, VectorType & v){
    typedef typename VectorType::ComponentType ScalarType;
    SymmetricMatrixType & D=M;
    VectorType & w = v;
    
    D = M.Inverse();
    const VectorType Dv = D*v;
    const ScalarType vDv = v.ScalarProduct(Dv);
    w = Dv/sqrt(1+vDv);
    D -= SymmetricMatrixType::RankOneTensor(w);
}

// Half disk metric (D,w), with in place adjustment of w
template<typename SymmetricMatrixType, typename VectorType, typename ScalarType>
void HalfDisk(ScalarType r, ScalarType delta, SymmetricMatrixType & D, VectorType & v){
    typedef SymmetricMatrixType Sym;
    D = (Sym::RankOneTensor(v)*(square(delta)-1)
         +Sym::Identity()*v.SquaredNorm())
    *square(r);
    v *= sqrt(std::max(0.,1.-square(delta*r)));
}
    
    /* // Older construction
     const MetricElementType & met = (*pDualMetric)(index);
     const VectorType & v = met.first/param.gridScale;
     const ScalarType & r = met.second;
     Voronoi1Vec<ReductionType>(&stencil.forward[0],
     v * sqrt(std::max(0.,1.-square(eps*r))),
     epsForward);
     typedef ReductionType::SymmetricMatrixType Sym;
     const Sym m =
     (Sym::RankOneTensor(v)*(square(eps)-1)
     +Sym::Identity()*v.SquaredNorm())
     *square(r);
     Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m);*/
    
}

// ------- Asymmetric quadratic model --------

template<size_t VDimension>
struct TraitsAsymmetricQuadratic : TraitsBase<VDimension> {
    typedef TraitsBase<VDimension> Superclass;
    Redeclare1Type(FromSuperclass,DiscreteType)
    Redeclare1Constant(FromSuperclass,Dimension)
    typedef typename Superclass::template Difference<0> DifferenceType;
    static const DiscreteType SymDimension = (Dimension*(Dimension+1))/2;
    static const DiscreteType nSymmetric=SymDimension, nForward=SymDimension;
};

template<size_t VDimension>
struct StencilAsymmetricQuadratic : HamiltonFastMarching<TraitsAsymmetricQuadratic<VDimension> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsAsymmetricQuadratic<VDimension> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare5Types(FromHFM,ParamDefault,ScalarType,Traits,VectorType,IndexType)
    Redeclare4Types(FromHFM,StencilType,ParamInterface,HFMI,IndexCRef)
    Redeclare1Constant(FromHFM,Dimension)
    ParamDefault param;
    ScalarType eps = 0.3, epsSide = 0.2;
    
    typedef typename Traits::template BasisReduction<Dimension> ReductionType;
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef LinearAlgebra::VectorPair<SymmetricMatrixType, VectorType> MetricElementType;
    typedef typename Traits::template DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pMetric;
    std::unique_ptr<MetricType> pDualMetric;

    typedef LinearAlgebra::VectorPair<VectorType, ScalarType> HalfDiskElementType;
    typedef typename Traits::template DataSource<HalfDiskElementType>  HalfDiskType;
    std::unique_ptr<HalfDiskType> pHalfDisk;


    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        MetricElementType met;
        if(pMetric){
            met = (*pMetric)(index);
            AsymmetricNorm::Dualize(met.first,met.second);
        } else if(pHalfDisk){
            const HalfDiskElementType & halfDisk = (*pHalfDisk)(index);
            met.second = halfDisk.first;
            AsymmetricNorm::HalfDisk(halfDisk.second,epsSide,met.first,met.second);
        } else {
            assert(pDualMetric);
            met = (*pDualMetric)(index);
        }
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],met.first);
        Voronoi1Vec<ReductionType>(&stencil.forward[0],met.second,eps);
        
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff: stencil.symmetric) diff.baseWeight*=hm2;
        for(auto & diff: stencil.forward)   diff.baseWeight*=hm2;
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        auto & io = that->io;
        Superclass::Setup(that);
        param.Setup(that);
        eps = io.template Get<ScalarType>("eps",eps);
        if(io.HasField("dualMetric")){
            pDualMetric = that->template GetField<MetricElementType>("dualMetric");
        } else if(io.HasField("metric")) {
            pMetric = that->template GetField<MetricElementType>("dualMetric");
        } else {
            pHalfDisk = that->template GetField<HalfDiskElementType>("halfDisk");
            epsSide=that->io.template Get<ScalarType>("epsForward",eps/1.5);
        }
    }
};


// ------- Three dimensional model ---------
/*
struct TraitsAsymmetricQuadratic3 : TraitsBase<3> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=6, nForward=6;
};

constexpr decltype(TraitsAsymmetricQuadratic3::boundaryConditions) TraitsAsymmetricQuadratic3::boundaryConditions;

struct StencilAsymmetricQuadratic3 : HamiltonFastMarching<TraitsAsymmetricQuadratic3>::StencilDataType {
    typedef HamiltonFastMarching<TraitsAsymmetricQuadratic3> HFM;
    typedef HFM::StencilDataType Superclass; // Copy pasted code below
    HFM::ParamDefault param;
    ScalarType eps = 0.3, epsSide = 0.2;
    
    typedef Traits::BasisReduction<Dimension> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef LinearAlgebra::VectorPair<SymmetricMatrixType, VectorType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pMetric;
    std::unique_ptr<MetricType> pDualMetric;
    
    typedef LinearAlgebra::VectorPair<VectorType, ScalarType> HalfDiskElementType;
    typedef Traits::DataSource<HalfDiskElementType>  HalfDiskType;
    std::unique_ptr<HalfDiskType> pHalfDisk;
    
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        MetricElementType met;
        if(pMetric){
            met = (*pMetric)(index);
            AsymmetricNorm::Dualize(met.first,met.second);
        } else if(pHalfDisk){
            const HalfDiskElementType & halfDisk = (*pHalfDisk)(index);
            met.second = halfDisk.first;
            AsymmetricNorm::HalfDisk(halfDisk.second,epsSide,met.first,met.second);
        } else {
            assert(pDualMetric);
            met = (*pDualMetric)(index);
        }
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],met.first);
        Voronoi1Vec<ReductionType>(&stencil.forward[0],met.second,eps);
        
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff: stencil.symmetric) diff.baseWeight*=hm2;
        for(auto & diff: stencil.forward)   diff.baseWeight*=hm2;
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        auto & io = that->io;
        Superclass::Setup(that);
        param.Setup(that);
        eps = io.Get<ScalarType>("eps",eps);
        if(io.HasField("dualMetric")){
            pDualMetric = that->GetField<MetricElementType>("dualMetric");
        } else if(io.HasField("metric")) {
            pMetric = that->GetField<MetricElementType>("dualMetric");
        } else {
            pHalfDisk = that->GetField<HalfDiskElementType>("halfDisk");
            epsSide=that->io.Get<ScalarType>("epsForward",eps/1.5);
        }
    }
};
*/
// ------- Three dimensional lifted model with one additional radius dimension ---------

struct TraitsAsymmetricQuadratic3p1 : TraitsBase<4> {
    typedef Difference<0> DifferenceType;
    constexpr static std::array<Boundary, Dimension> boundaryConditions =
    {{Boundary::Closed, Boundary::Closed, Boundary::Closed, Boundary::Closed}};
    static const DiscreteType nSymmetric=6+1, nForward=6;
    
    // Stencils actually depend on all coordinates (no multiplier). This is to get proper domain parametrization.
    static const DiscreteType nStencilDependencies=1;
    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{3}};
};

constexpr decltype(TraitsAsymmetricQuadratic3p1::boundaryConditions)  TraitsAsymmetricQuadratic3p1::boundaryConditions;
constexpr decltype(TraitsAsymmetricQuadratic3p1::stencilDependencies) TraitsAsymmetricQuadratic3p1::stencilDependencies;

struct StencilAsymmetricQuadratic3p1 : HamiltonFastMarching<TraitsAsymmetricQuadratic3p1>::StencilDataType {
    typedef HamiltonFastMarching<TraitsAsymmetricQuadratic3p1> HFM;
    typedef HFM::StencilDataType Superclass; // Copy pasted code below
    HFM::ParamDefault param;
    ScalarType eps = 0.3, epsSide = 0.2;
    
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::SymmetricMatrixType SymmetricMatrixType;
    typedef ReductionType::VectorType VectorType; // Redefinition of VectorType as a 3d vector
    typedef LinearAlgebra::VectorPair<SymmetricMatrixType, VectorType> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pMetric;
    std::unique_ptr<MetricType> pDualMetric;
    
    typedef LinearAlgebra::VectorPair<VectorType, ScalarType> HalfDiskElementType;
    typedef Traits::DataSource<HalfDiskElementType>  HalfDiskType;
    std::unique_ptr<HalfDiskType> pHalfDisk;
    
    typedef Traits::DataSource<ScalarType> BundleCostType;
    std::unique_ptr<BundleCostType> pBundleCost;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        MetricElementType met;
        if(pMetric){
            met = (*pMetric)(index);
            AsymmetricNorm::Dualize(met.first,met.second);
        } else if(pHalfDisk){
            const HalfDiskElementType & halfDisk = (*pHalfDisk)(index);
            met.second = halfDisk.first;
            AsymmetricNorm::HalfDisk(halfDisk.second,epsSide,met.first,met.second);
        } else {
            assert(pDualMetric);
            met = (*pDualMetric)(index);
        }
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0],met.first);
        Voronoi1Vec<ReductionType>(&stencil.forward[0],met.second,eps);
        
        const ScalarType hm2 = 1/square(param.gridScale);
        for(auto & diff: stencil.symmetric) diff.baseWeight*=hm2;
        for(auto & diff: stencil.forward)   diff.baseWeight*=hm2;
        
        stencil.symmetric[6].offset = OffsetType::Constant(0);
        stencil.symmetric[6].offset[3]=1;
        stencil.symmetric[6].baseWeight=square((*pBundleCost)(index)/param.dependScale);
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        auto & io = that->io;
        Superclass::Setup(that);
        param.Setup(io,io.Get<ScalarType>("bundleScale",1));
        eps = io.Get<ScalarType>("eps",eps);
        if(io.HasField("dualMetric")){
            pDualMetric = that->GetField<MetricElementType>("dualMetric");
        } else if(io.HasField("metric")) {
            pMetric = that->GetField<MetricElementType>("dualMetric");
        } else {
            pHalfDisk = that->GetField<HalfDiskElementType>("halfDisk");
            epsSide=that->io.Get<ScalarType>("epsForward",eps/1.5);
        }
        pBundleCost = that->GetField<ScalarType>("bundleCost");
    }
};

/*
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps = 0.2, epsForward = 0.3;
    typedef Traits::BasisReduction<3> ReductionType;
    typedef ReductionType::VectorType Vector3;
    /// Physical speed, bundle speed.
    typedef LinearAlgebra::Vector<ScalarType, 2> SpeedPair;
    typedef LinearAlgebra::VectorPair<Vector3, SpeedPair> MetricElementType;
    typedef Traits::DataSource<MetricElementType> MetricType;
    std::unique_ptr<MetricType> pDualMetric;
    
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        const MetricElementType & met = (*pDualMetric)(index);
        const Vector3 & v = met.first/param.gridScale;
        const ScalarType & r = met.second[0];
        
        Voronoi1Vec<ReductionType>(&stencil.forward[0],
                                   v * sqrt(std::max(0.,1.-square(eps*r))),
                                   epsForward);
        typedef ReductionType::SymmetricMatrixType Sym;
        const Sym m =
        (Sym::RankOneTensor(v)*(square(eps)-1)
         +Sym::Identity()*v.SquaredNorm())
        *square(r);
        Voronoi1Mat<ReductionType>(&stencil.symmetric[0], m);
        
        for(int i=0; i<6; ++i){
            stencil.symmetric[i].offset[3]=0;
            stencil.forward[i].offset[3]=0;}
        
        stencil.symmetric[6].offset = OffsetType::Constant(0);
        stencil.symmetric[6].offset[3]=1;
        stencil.symmetric[6].baseWeight=square(met.second[1]/param.gridScale);
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        param.Setup(that);
        eps=that->io.Get<ScalarType>("eps",eps);
        epsForward=that->io.Get<ScalarType>("epsForward",eps*1.5);
        pDualMetric = that->GetField<MetricElementType>("dualMetric");
    }
};*/


#endif /* AsymmetricQuadratic_h */
