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
    typedef int DiscreteType;
    typedef int_least8_t ShortType;
    static constexpr ScalarType Infinity() {return std::numeric_limits<ScalarType>::infinity();}
    static constexpr ScalarType mathPi = M_PI;
    
    typedef LinearAlgebra::Point<ScalarType, Dimension> PointType;
    typedef LinearAlgebra::Vector<ScalarType, Dimension> VectorType;
    typedef LinearAlgebra::Point<DiscreteType, Dimension> IndexType;
    typedef LinearAlgebra::Vector<ShortType, Dimension> OffsetType;
    template<size_t n> using
    BasisReduction=LinearAlgebra::BasisReduction<ScalarType, DiscreteType, n>;
    
    template<typename T, size_t n> using Array=LinearAlgebra::Array<T,n>;
    
    // The following stencil specs should be overriden in subclass
    static const DiscreteType
    nForward=0, nSymmetric=0,
    nMax=1, nMaxForward=0, nMaxSymmetric=0;
    
    // How many scalars determine the weights ?
    template<size_t VMultSize, typename Dummy=void> struct Difference;
    
    // Data sources
    template<typename E> struct DataSource {
        virtual ~DataSource(){};
//        typedef typename std::conditional<std::numeric_limits<E>::is_specialized, E, const E &>::type ReturnType;
        typedef E ReturnType;
        virtual bool CheckDims(const IndexType &) const = 0;
        virtual ReturnType operator()(const IndexType &) const = 0;
    };
};

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

/* // GCC cannot match the following template in a template, so a macro is used instead
template<int VDim,int k> std::ostream & operator <<
(std::ostream & os, const typename TraitsBase<VDim>::template Difference<k> & diff) {
     return os << "{" << diff.baseWeight << "," << diff.offset  << "}";}*/

#define __PrintDiff(__VDim) \
template<size_t k> std::ostream & operator << \
(std::ostream & os, const TraitsBase<__VDim>::Difference<k> & diff) { \
    return os << "{" << diff.baseWeight << "," << diff.offset  << "}";}

__PrintDiff(2) __PrintDiff(3) __PrintDiff(4) __PrintDiff(5)

// --------------------- Tensor decomposition ---------------------

template<typename ReductionType, int VDimShift=0, typename DiffType,
typename SymmetricMatrixType = typename ReductionType::SymmetricMatrixType,
typename ScalarType = typename ReductionType::ScalarType>
DiffType * Voronoi1Mat(DiffType * pDiff,
                       const SymmetricMatrixType & d,
                       ScalarType tol=1000*std::numeric_limits<ScalarType>::epsilon()){
    // Get the tensor decomposition.
    const auto & decomp = ReductionType::TensorDecomposition(d);
    // Fill in the scheme, omitting excessively low values
    for(const auto & offsetWeight : decomp){
        const auto & offset = offsetWeight.first;
        const ScalarType & weight = offsetWeight.second;
        pDiff->baseWeight = weight<tol ? 0. : weight;
        pDiff->offset.fill(0);
        for(int i=0; i<offset.size(); ++i)
            pDiff->offset[VDimShift+i]=offset[i];
        ++pDiff;
    }
    return pDiff;
}

template<typename ReductionType, int VDimShift=0, typename DiffType,
typename VectorType = typename ReductionType::VectorType,
typename ScalarType = typename ReductionType::ScalarType
>
DiffType * Voronoi1Vec(DiffType * pDiff,
                       const VectorType & v, ScalarType eps,
                       ScalarType tol=1000*std::numeric_limits<ScalarType>::epsilon()){
    // Build tensor
    typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
    const SymmetricMatrixType d =
    SymmetricMatrixType::RankOneTensor(v)*(1.-square(eps))
    +v.SquaredNorm()*SymmetricMatrixType::Identity()*square(eps);
    // Get decomposition
    const auto & decomp = ReductionType::TensorDecomposition(d);
    // Fill in the scheme, omitting excessively low values, with reorientation
    for(const auto & offsetWeight : decomp){
        const auto & offset = offsetWeight.first;
        const ScalarType & weight = offsetWeight.second;
        
        const ScalarType scal = v.ScalarProduct(VectorType::CastCoordinates(offset));
        pDiff->baseWeight = (weight<tol || std::abs(scal)<tol) ? 0. : weight;
        
        pDiff->offset.fill(0);
        for(int i=0; i<offset.size(); ++i) {
            pDiff->offset[VDimShift+i]= scal>=0 ? offset[i] : -offset[i];}
        ++pDiff;
    }
    return pDiff;
};


#endif /* CommonTraits_h */
