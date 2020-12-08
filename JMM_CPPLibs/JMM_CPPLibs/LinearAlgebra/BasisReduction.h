// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

/*
 Notice to users.
 If you use these methods in an academic or commercial project, it would be kind to cite at least one of the following papers, where the applications of basis reduction to PDE discretization were first introduced.
 - (On anisotropic fast marching)
J.-M. Mirebeau, “Anisotropic Fast-Marching on cartesian grids using Lattice Basis Reduction,” SIAM J. Numer. Anal., vol. 52, no. 4, pp. 1573–1599, Jan. 2014.
 - (On anisotropic diffusion)
J. Fehrenbach and J.-M. Mirebeau, “Sparse Non-negative Stencils for Anisotropic Diffusion,” Journal of Mathematical Imaging and Vision, pp. 1–25, 2013.
 */

#ifndef BasisReduction_h
#define BasisReduction_h

#include "SymmetricMatrixType.h"
#include "../Macros/RedeclareTypes.h"
#include "../Macros/DependentFalse.h"
#include "../Macros/ExportArrow.h"

namespace LinearAlgebra {
    
template<typename TScalar, typename TDiscrete, size_t VDimension>
struct BasisReduction {
    typedef TScalar ScalarType;
    typedef TDiscrete DiscreteType;
    static const size_t Dimension = VDimension;
    typedef Vector<ScalarType, Dimension> VectorType;
    typedef Vector<DiscreteType, Dimension> DiscreteVectorType;
    typedef SymmetricMatrix<ScalarType, Dimension> SymmetricMatrixType;
    
    typedef std::array<DiscreteVectorType, Dimension>   BasisType;
    typedef std::array<DiscreteVectorType, Dimension+1> SuperbaseType;
    static const size_t SymDimension = (Dimension*(Dimension+1))/2;
	static constexpr int KKTDimension = SymDimension; // Differs in dimension 4, see VoronoiReduction
    struct TensorDecompositionType {
        std::array<DiscreteVectorType,KKTDimension> offsets;
        std::array<ScalarType,KKTDimension> weights;
    };
    
    static const BasisType CanonicalBasis();
    static const SuperbaseType CanonicalSuperBase();
    static void ReducedBasis(const SymmetricMatrixType &, BasisType &); // Minkowksi reduced basis (Dimension=2 only)
    static void ObtuseSuperbase(const SymmetricMatrixType &, SuperbaseType &); // Selling reduced basis (Dimension<=3)
    static TensorDecompositionType TensorDecomposition(const SymmetricMatrixType &);
    static int maxIt;
	
	/// Complementary indices for set {0,1,..., Dimension}
	static std::array<int,Dimension-1> ComplementIndices(int,int);
	/// Linearized index for pair (i,j), 0<=i<j<=Dimension, fast j
	static int LinearizeIndices(int,int);
	
	/// Structure dedicated to enumerating the obtuse superbases along the line between two symmetric matrices
	struct SellingPath {
		SymmetricMatrixType D0,D1;
		SuperbaseType sb;
		std::array<DiscreteVectorType,SymDimension> offsets;
		std::array<ScalarType,SymDimension> weights0,weights1;
		SellingPath(const SymmetricMatrixType &, const SymmetricMatrixType &, ScalarType);
		
		// Next step, with a lower and upper bound.
		// Returns pair corresponding to the previously obtuse pair.
		std::pair<int,int> NextStep(ScalarType &) const;
		void MakeStep(const std::pair<int,int> &);
		PrintSelfMacro(SellingPath);
	};
protected:
    template<size_t VD, typename Dummy> struct TensorDecompositionHelper_;
	typedef TensorDecompositionHelper_<Dimension, void> TensorDecompositionHelper;
};
	
template<typename TS,typename TD, size_t VD> int BasisReduction<TS,TD,VD>::maxIt=200;

#include "Implementation/BasisReduction.hpp"
}


#endif /* BasisReduction_h */
