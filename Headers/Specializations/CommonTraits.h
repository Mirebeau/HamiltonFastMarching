// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef CommonTraits_h
#define CommonTraits_h

#include "Base/HFMInterface.h"

#include "JMM_CPPLibs/LinearAlgebra/BasisReduction.h"
#include "JMM_CPPLibs/LinearAlgebra/SquareCube.h"

#include "Base/BaseGrid.h"
#include "Base/PeriodicGrid.h" //Boundary_AllClosed.
#include "Base/EulerianStencil.h"

#ifdef HighVoronoi // Will need high dimensional Voronoi reduction
#include <type_traits> // std::conditional
#include "JMM_CPPLibs/LinearAlgebra/VoronoiReduction.h"
#endif

template<int VDim> struct TraitsBase {
    static constexpr int Dimension = VDim;
    typedef BaseGrid<VDim, double, int, int_least8_t> BaseDomain;
    Redeclare6Types(BaseDomain,DiscreteType,IndexType,IndexDiff,IndexCRef,ShortType,OffsetType)
    Redeclare3Types(BaseDomain,ScalarType,PointType,VectorType)
    template<typename T, size_t n> using Array=typename BaseDomain::template Array<T,n>;
    
    static constexpr ScalarType Infinity() {return std::numeric_limits<ScalarType>::infinity();}
    static constexpr ScalarType mathPi =  3.141592653589793238462643383279502884L;
    
    // Tensor decomposition based on Voronoi first reduction.
    // (Involved in many PDE discretization schemes.)
#ifdef HighVoronoi // Change the tensor decomposition in high dimension
    template<size_t n> using BasisReduction=typename std::conditional<
    n<=3,
    LinearAlgebra::BasisReduction<ScalarType, DiscreteType, n>,
    VoronoiFirstReduction<ScalarType,Dimension>
    >::type;
#else
    template<size_t n> using
    BasisReduction=LinearAlgebra::BasisReduction<ScalarType, DiscreteType, n>;
#endif
        
    /** On which coordinates do the adaptive stencils depend (default : none)
     Only applies if there is a multiplier field*/
    static const DiscreteType nStencilDependencies=0;
    typedef std::array<DiscreteType, nStencilDependencies> StencilDepType_None;
    constexpr static const StencilDepType_None stencilDependencies = {{}};
    
    /// Boundary conditions for the PDE (default : closed box)
    constexpr static const Boundary_AllClosed boundaryConditions{};
    
    /// A DataSource is a pure interface, multi-dimensional-array-like for providing pointwise data on a grid
    template<typename E> struct DataSource {
        virtual ~DataSource(){};
//        typedef typename std::conditional<std::numeric_limits<E>::is_specialized, E, const E &>::type ReturnType;
        typedef E ReturnType;
        virtual bool CheckDims(const IndexType &) const = 0;
        virtual ReturnType operator()(const IndexType &) const = 0;
    };
    
	// Likely hidden in child class
	using DistanceGuess = LinearAlgebra::SymmetricMatrix<ScalarType,Dimension>;
};
// Linker wants the following two lines for some obscure reason.
template<int VD> constexpr const typename TraitsBase<VD>::StencilDepType_None TraitsBase<VD>::stencilDependencies;
template<int VD> constexpr const Boundary_AllClosed TraitsBase<VD>::boundaryConditions;


// --------------------- Tensor decomposition ---------------------

/** PDE discretization helper.
 This function leverages Voronoi's first reduction of quadratic forms to decompose
 a symmetric positive definite matrix into a positively weighted sum of
 rank one tensors with integer coefficients. D = sum_i lambda_i e_i x e_i.
 The (lambda_i, e_i) are returned as baseweight and offset of differences (pDiff).
 Excessively small lambda_i are erased, based on tol parameter.
 See : Jean-Marie Mirebeau,
 Anisotropic fast-marching on cartesian grids using Voronoi’s first reduction of quadratic forms,
 (preprint available on Arxiv) */
template<typename ReductionType, int VDimShift=0, typename DiffType,
typename SymmetricMatrixType = typename ReductionType::SymmetricMatrixType,
typename ScalarType = typename ReductionType::ScalarType>
DiffType * Voronoi1Mat(DiffType * pDiff,
                       const SymmetricMatrixType & d,
                       ScalarType tol=1000*std::numeric_limits<ScalarType>::epsilon()){
    const ScalarType atol = tol*d.Trace(); // Absolute tolerance
    
    // Get the tensor decomposition.
    const auto & decomp = ReductionType::TensorDecomposition(d);
    // Fill in the scheme, omitting excessively low values
    auto offsetIt = decomp.offsets.begin();
    auto weightIt = decomp.weights.begin();
    for(; weightIt!=decomp.weights.end(); ++weightIt, ++offsetIt){
        const ScalarType weight = *weightIt;
        const auto & offset = *offsetIt;
        pDiff->baseWeight = weight<atol ? 0. : weight;
        pDiff->offset.fill(0);
        for(int i=0; i< offset.size(); ++i)
            pDiff->offset[VDimShift+i]= offset[i];
        ++pDiff;
    }
    return pDiff;
}

/** PDE discretization helper.
 This function leverages Voronoi's first reduction of quadratic forms to
 construct approximations <g,v>_+^2 = sum_i lambda_i <g, -e_i>_+^2,
 where _+ denotes positive part, lambda_i>=0, e_i is an offset.
 
 This is used in the finite differences approximation
 <grad u(x),v>_+^2 ~ sum_i lambda_i (g(x) - g(x+e_i))_+^2
 
 See : Jean-Marie Mirebeau,
 Fast Marching methods for Curvature Penalized Shortest Paths, 2017
*/
template<typename ReductionType, int VDimShift=0> struct Voronoi1Vec {
    typedef typename ReductionType::VectorType VectorType;
    typedef typename ReductionType::ScalarType ScalarType;
    static const int Dimension = VectorType::Dimension;
    ScalarType
    eps=0.1, //
    angleTol=1./sqrt(ScalarType(Dimension-1)), // Eliminate offsets e s.t. |e| |v| > sqrt(1+angleTol^2) <e,v>
    weightTol=1000*std::numeric_limits<ScalarType>::epsilon(); // Eliminate excessively small weights.
    
    template<typename DiffType>
    DiffType * operator()(DiffType * pDiff, const VectorType & v){
        // Build tensor
        typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
        const SymmetricMatrixType exact = SymmetricMatrixType::RankOneTensor(v);
        const SymmetricMatrixType relaxed =
        exact*(1.-square(eps))
        +v.SquaredNorm()*SymmetricMatrixType::Identity()*square(eps);
        // Get decomposition
        const auto & decomp = ReductionType::TensorDecomposition(relaxed);
        
        auto offsetIt = decomp.offsets.begin();
        auto weightIt = decomp.weights.begin();
        // Fill in the scheme, omitting excessively low values, with reorientation
        const ScalarType vNorm2 = v.SquaredNorm();
        const ScalarType wTol = vNorm2*weightTol; // Absolute tolerance on weights
        const ScalarType angTol = (1+square(angleTol))/vNorm2;

        for(; weightIt!=decomp.weights.end(); ++weightIt, ++offsetIt){
            const ScalarType weight = *weightIt;
            const auto & offset = *offsetIt;
            
            const VectorType off = VectorType::CastCoordinates(offset);
            
            const ScalarType scal = v.ScalarProduct(off);
            pDiff->baseWeight = (weight<wTol || square(scal)*angTol<off.SquaredNorm()) ? 0. : weight;
            
            pDiff->offset.fill(0);
            for(int i=0; i<offset.size(); ++i) {
                pDiff->offset[VDimShift+i]= scal<=0 ? offset[i] : -offset[i];} // Note : <v,e_i> <= 0
            ++pDiff;
        }
        return pDiff;
    }
};

#endif /* CommonTraits_h */
