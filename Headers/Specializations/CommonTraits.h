// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef CommonTraits_h
#define CommonTraits_h


#include "LinearAlgebra/SymmetricMatrixType.h"
#include "LinearAlgebra/ArrayType.h"
#include "LinearAlgebra/BasisReduction.h"
#include "LinearAlgebra/SquareCube.h"

template<int VDim> struct TraitsBase {
    static const int Dimension = VDim;
    typedef double ScalarType;
    typedef int DiscreteType; /// Used for array indexing and multi-index components.
    typedef int_least8_t ShortType;
    static constexpr ScalarType Infinity() {return std::numeric_limits<ScalarType>::infinity();}
    static constexpr ScalarType mathPi =  3.141592653589793238462643383279502884L;
    
    typedef LinearAlgebra::Point<ScalarType, Dimension> PointType;
    typedef LinearAlgebra::Vector<ScalarType, Dimension> VectorType;
    typedef LinearAlgebra::Point<DiscreteType, Dimension> IndexType; /// Multi-index of a grid point
    typedef LinearAlgebra::Vector<ShortType, Dimension> OffsetType; /// Used for offsets between neighbor indices
    template<size_t n> using
    BasisReduction=LinearAlgebra::BasisReduction<ScalarType, DiscreteType, n>; // Ingredient of many PDE discretization schemes.
    
    template<typename T, size_t n> using Array=LinearAlgebra::Array<T,n>;
    
    /** A Difference is a basic component of a PDE scheme. It is the data of an offset and weight.
     The weight which is either specified directly or as a baseweight and a multiplier index, within [0,VMultSize[.*/
    template<size_t VMultSize, typename Dummy=void> struct Difference;

    
    /** A PDE discretization is specified via a set of differences, of several types, 
     which number is fixed in advance by the following parameters.
     These (empty) defaults need to be redefined (=shadowed) in a subclass.*/
    static const DiscreteType
    nForward=0, nSymmetric=0,
    nMax=1, nMaxForward=0, nMaxSymmetric=0;
    
    /// On which coordinates do the adaptive stencils depend (default : none)
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
};

template<int VD> constexpr const typename TraitsBase<VD>::StencilDepType_None TraitsBase<VD>::stencilDependencies;
template<int VD> constexpr const Boundary_AllClosed TraitsBase<VD>::boundaryConditions;

template<int VDim> template<size_t VMultSize, typename Dummy> struct TraitsBase<VDim>::Difference {
    typedef TraitsBase<VDim> T; typedef T::ScalarType ScalarType; typedef T::OffsetType OffsetType;
    typedef T::ShortType ShortType;
        ScalarType baseWeight; OffsetType offset; ShortType multIndex; static constexpr size_t multSize=VMultSize;
        typedef LinearAlgebra::Point<ScalarType, multSize> MultiplierType;
        ScalarType Weight(const MultiplierType & mult) const {return baseWeight*square(mult[multIndex]);}
//        PrintSelfMacro(Difference);
};

template<int VDim> template<typename Dummy> struct TraitsBase<VDim>::Difference<1,Dummy> {
    typedef TraitsBase<VDim> T; typedef T::ScalarType ScalarType; typedef T::OffsetType OffsetType;
    ScalarType baseWeight; OffsetType offset; static constexpr size_t multSize=1;
    typedef ScalarType MultiplierType;
    ScalarType Weight(const MultiplierType & mult) const {
        return baseWeight*square(mult);}
};

template<int VDim> template<typename Dummy> struct TraitsBase<VDim>::Difference<0,Dummy> {
    typedef TraitsBase<VDim> T; typedef T::ScalarType ScalarType; typedef T::OffsetType OffsetType;
    ScalarType baseWeight; OffsetType offset; static constexpr size_t multSize=0;
    struct MultiplierType {};
    ScalarType Weight(MultiplierType) const {return baseWeight;}
};

/* // Print differences
 // GCC cannot match the following template in a template, so a macro is used instead
template<int VDim,int k> std::ostream & operator <<
(std::ostream & os, const typename TraitsBase<VDim>::template Difference<k> & diff) {
     return os << "{" << diff.baseWeight << "," << diff.offset  << "}";}*/

#define __PrintDiff(__VDim) \
template<size_t k> std::ostream & operator << \
(std::ostream & os, const TraitsBase<__VDim>::Difference<k> & diff) { \
    return os << "{" << diff.baseWeight << "," << diff.offset  << "}";}

__PrintDiff(2) __PrintDiff(3) __PrintDiff(4) __PrintDiff(5)

// --------------------- Tensor decomposition ---------------------

/** PDE discretization helper.
 This function leverages Voronoi's first reduction of quadratic forms to decompose
 a symmetric positive definite matrix into a positively weighted sum of
 rank one tensors with integer coefficients. D = sum_i lambda_i e_i x e_i.
 The (lambda_i, e_i) are returned as baseweight and offset of differences (pDiff).
 Excessively small lambda_i are erased, based on tol parameter.
 See : Jean-Marie Mirebeau, (preprint available on Arxiv)
 Anisotropic fast-marching on cartesian grids using Voronoi’s first reduction of quadratic forms */
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
 Fast Marching methods for Curvature Penalized Shortest Paths,
 (preprint available on Arxiv) */
template<typename ReductionType, int VDimShift=0, typename DiffType,
typename VectorType = typename ReductionType::VectorType,
typename ScalarType = typename ReductionType::ScalarType
>
DiffType * Voronoi1Vec(DiffType * pDiff,
                       const VectorType & v, ScalarType eps,
                       ScalarType tol=1000*std::numeric_limits<ScalarType>::epsilon()){
    const ScalarType atol = v.SquaredNorm()*tol; // Absolute tolerance
    
    // Build tensor
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    const SymmetricMatrixType d =
    SymmetricMatrixType::RankOneTensor(v)*(1.-square(eps))
    +v.SquaredNorm()*SymmetricMatrixType::Identity()*square(eps);
    // Get decomposition
    const auto & decomp = ReductionType::TensorDecomposition(d);
    auto offsetIt = decomp.offsets.begin();
    auto weightIt = decomp.weights.begin();
    // Fill in the scheme, omitting excessively low values, with reorientation
    for(; weightIt!=decomp.weights.end(); ++weightIt, ++offsetIt){
        const ScalarType weight = *weightIt;
        const auto & offset = *offsetIt;

        const ScalarType scal = v.ScalarProduct(VectorType::CastCoordinates(offset));
        pDiff->baseWeight = (weight<atol || std::abs(scal)<atol) ? 0. : weight;
        
        pDiff->offset.fill(0);
        for(int i=0; i<offset.size(); ++i) {
            pDiff->offset[VDimShift+i]= scal<=0 ? offset[i] : -offset[i];} // Note : <v,e_i> <= 0
        ++pDiff;
    }
    return pDiff;
};


#endif /* CommonTraits_h */
