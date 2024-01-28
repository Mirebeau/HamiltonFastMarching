//
//  VoronoiReduction.h
//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 06/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef VoronoiReduction_h
#define VoronoiReduction_h

/*
 The purpose of this program is to implement Voronoi's first reduction in dimension 4 and 5.
 The steps are :
 - (optional) A basis reduction step, similar to the LLL algorithm.
 - Finding the optimum in Voronoi's linear problem, using a simplex like algorithm.
 Note that some of the vertices has high multiplicity (64 and 400 respectively).
 - Computing the Karush-Kuhn-Tucker relations, which requires solving a low dimensional linear program at the vertices with high multiplicity.
 
 In the first two steps, one must keep track of the change of coordinates.

 This code contains a significant amount of raw data, preprocessed using the Mathematica(r) software, and related to the theory of "perfect quadratic forms".
 */


#include "../Macros/ExportArrow.h"
#include "SymmetricMatrixType.h"

#ifdef Voronoi6
// Various macros and typedefs to allow executing the "CUDA" code on the CPU.
#define Scalar_macro
#define DOUBLE
typedef double Scalar;
#define __restrict__
#include <cmath>
#include <iostream>
using std::sqrt; using std::max; using std::min;
#define IOSTREAM
#define SIMPLEX_TOL 2e-14

/* The 6-dimensional Voronoi reduction algorithm is implemented in a different library, see
https://github.com/Mirebeau/AdaptiveGridDiscretizations
If you do not need it, then please #undef Voronoi6
If you need it, the please put the correct include path below
*/
#ifdef _MSC_VER 
#include "C:/Users/jmmir/Documents/GitHub/AdaptiveGridDiscretizations/agd/Eikonal/HFM_CUDA/cuda/Geometry6.h" // My Dell
#else
#include "/Users/jean-mariemirebeau/Dropbox/Programmes/GithubM1/AdaptiveGridDiscretizations/agd/Eikonal/HFM_CUDA/cuda/Geometry6.h" // My Macbook M1
//#include "Eikonal/HFM_CUDA/CUDA/geometry6.h"
//#include "/Users/mirebeau/Dropbox/Programmes/Github/AdaptiveGridDiscretizations/agd/Eikonal/HFM_CUDA/cuda/Geometry6.h"
#endif

#endif

#include "Implementation/LinProg/Siedel_Hohmeyer_LinProg.h"



template<typename TScalar,int VDimension>
struct VoronoiFirstReduction {
    typedef TScalar ScalarType;
    static const int Dimension=VDimension;
    typedef LinearAlgebra::Vector<ScalarType,Dimension> VectorType;
    typedef LinearAlgebra::Matrix<ScalarType,Dimension,Dimension> MatrixType;
    typedef LinearAlgebra::SymmetricMatrix<ScalarType,Dimension> SymmetricMatrixType;
    static constexpr ScalarType infinity = std::numeric_limits<ScalarType>::infinity();
//    static constexpr int badIndex = std::numeric_limits<int>::max(); // Not constexpr ?
    static constexpr int badIndex = 1234567890;
    
    struct SimplexStateType {
        SymmetricMatrixType m;
        MatrixType a;
        int vertex;
        ScalarType objective;
        SimplexStateType(const SymmetricMatrixType & m0)
        :m(m0),a(MatrixType::Identity()),vertex(-1),objective(infinity){};
        PrintSelfMacro(SimplexStateType)
    };
    
    // Basis reduction. Optional argument is a tolerance parameter
    static void GreedyBasis(SimplexStateType &, ScalarType=0);
    // The basis is Minkowski reduced if the result is non-positive (not always the case)
    static ScalarType IsMinkowskiReduced(const SymmetricMatrixType &);
    
    // Simplex algorithm
    static void Minimize(SimplexStateType &);

    // Kuhn-Tucker relations
    typedef int DiscreteType;
    typedef LinearAlgebra::Vector<DiscreteType, Dimension> OffsetType;
    static constexpr int SymDimension = Dimension*(Dimension+1)/2;
    static constexpr int MatDimension = Dimension*Dimension;
	 // Select a special representative among all possible decompositions in dims 4 and 5
	static const bool useLipschitzDecomposition4 = true;
	static const bool useLipschitzDecomposition5 = true;
	static constexpr int KKTDimension = SymDimension + 2*int(Dimension==4 && useLipschitzDecomposition4);
    struct KKTRelationType {
        std::array<OffsetType,KKTDimension> offsets;
        std::array<ScalarType,KKTDimension> weights;
        PrintSelfMacro(KKTRelationType)
    };
    static KKTRelationType KKT(const SimplexStateType &);
    static KKTRelationType TensorDecomposition(const SymmetricMatrixType &,ScalarType=0);
    
protected:
    // Basis reduction
    static void SortBasis(SimplexStateType &);
    template<int> static int GreedyBasisStep(SimplexStateType & state, ScalarType);

    // Simplex algorithm
    static void FirstGuess(SimplexStateType &);
    static bool BetterNeighbor(SimplexStateType &);
    static ScalarType Scal(const SymmetricMatrixType &, const SymmetricMatrixType &);
    
    static const int nVertices = (Dimension<=3 ? 1 : Dimension==4 ? 2 : 3);
    static const SymmetricMatrixType & GetVertex(int);
    static constexpr int NNeighbors(int);
    static void SetNeighbor(SimplexStateType &, int);
};


#include "Implementation/VoronoiReduction_Basis.hxx"
#include "Implementation/VoronoiReduction_Minimize_Data.hxx"
#include "Implementation/VoronoiReduction_Minimize.hxx"
#include "Implementation/VoronoiReduction_KKT.hxx"

#endif /* VoronoiReduction_h */
