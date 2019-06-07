//
//  VoronoiReduction_Basis.hxx
//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 06/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef VoronoiReduction_Basis_hxx
#define VoronoiReduction_Basis_hxx

template<typename TS, int VD> void VoronoiFirstReduction<TS,VD>::
GreedyBasis(SimplexStateType & state, ScalarType tol) {
    SortBasis(state);
    // Loop until basis is reduced
    for(int d=1; d<Dimension;){
        if(Dimension>=2 && d==1) {d=GreedyBasisStep<1>(state,tol);}
        if(Dimension>=3 && d==2) {d=GreedyBasisStep<2>(state,tol);}
        if(Dimension>=4 && d==3) {d=GreedyBasisStep<3>(state,tol);}
        if(Dimension>=5 && d==4) {d=GreedyBasisStep<4>(state,tol);}
    }
}

template<typename TS, int VD> void VoronoiFirstReduction<TS,VD>::
SortBasis(SimplexStateType & state){
    const SymmetricMatrixType m = state.m;
    
    // Find ordering of norms
    std::array<int, Dimension> order;
    for(int d=0; d<Dimension; ++d) order[d]=d;
    std::sort(order.begin(), order.end(),
              [&m](int i, int j){return m(i,i)<=m(j,j);});
    
    // Apply to symmetric matrix
    for(int i=0; i<Dimension; ++i){
        for(int j=0; j<=i; ++j){
            state.m(i,j) = m(order[i],order[j]);}
    }
    
    // Apply to basis
    const MatrixType a = state.a;
    for(int i=0; i<Dimension; ++i){ // TODO : Check
        const int j = order[i];
        VectorType & rowi = *(VectorType*)&state.a(i,0);
        const VectorType & rowj = *(const VectorType*)&a(j,0);
        rowi = rowj;

    }
}

template<typename TS,int VD> template<int d> int VoronoiFirstReduction<TS,VD>::
GreedyBasisStep(SimplexStateType & state, ScalarType tol) {
    assert(d<Dimension);
    SymmetricMatrixType & m = state.m;
    MatrixType & a = state.a;
    
    // Subtract to the d-th row a good integral combination of the former.
    typedef LinearAlgebra::SymmetricMatrix<ScalarType,d> SymD;
    typedef LinearAlgebra::Vector<ScalarType,d> VecD;
    const SymD & md = *(const SymD*)&m;
    const VecD & vd = *(const VecD*)&m(d,0);
    
    // TODO : proper solve function
    VecD mInvV = md.CGSolve(vd); // md.Inverse()*vd;
    
    // Create the coefficients of the new basis vector
    // Remark, rounding is not enough to get a Minkowski basis at all times.
    // A bit of exhaustive search would be needed in principle.
    // TODO : try the two values, if rounding an entry close to 0.5.
    VectorType coef;
    bool untouched = true;
    for(int i=0; i<d; ++i) {
        const ScalarType c = -round(mInvV[i]);
        untouched &= (c==0);
        coef[i] = c;
    }
    
    
    /*{
        std::cout
        ExportVarArrow(mInvV) << "\n"
        ExportVarArrow(coef) << "\n"
        ExportVarArrow(m) << "\n"
        ;
    }*/
    
    if(untouched) return d+1; // Go to next step
    coef[d]=1;
    for(int i=d+1; i<Dimension; ++i) coef[i]=0;
    
    const ScalarType prevNorm2 = m(d,d);
    const VectorType mCoef = m*coef;
    const ScalarType nextNorm2 = mCoef.ScalarProduct(coef); //md1.SquaredNorm(coef);
    
    // If norm reduction is not substantial, then go to next step
    if(nextNorm2*(1+tol) >= prevNorm2) return d+1;
    
    
    /*{
        std::cout
        ExportVarArrow(m.Gram(a.Inverse().Transpose()))
        << "\n"
        ExportVarArrow(mInvV) << "\n"
        ExportVarArrow(coef) << "\n"
        ExportVarArrow(m) << "\n"
        ExportVarArrow(a) << "\n"
        ExportVarArrow(mCoef)
        << "\n";
//        m0 = m;
 //       std::cout ExportVarArrow
    }*/
    
    // Otherwise apply changes.
    // --> to the symmetric matrix
    for(int i=0; i<d; ++i) m(d,i) = mCoef[i];
    m(d,d) = nextNorm2;
    for(int i=d+1; i<Dimension; ++i) m(d,i) = mCoef[i];
    
    // --> to the basis
    // Due to matrix layout (row major), rows of "a" are regarded as vectors
    VectorType & rowd = *(VectorType*)&a(d,0);
    for(int i=0; i<d; ++i) {
        VectorType & rowi = *(VectorType*)&a(i,0);
        rowd += coef[i] * rowi;
    }

    /*
    {
        std::cout << "Changes applied\n"
        ExportVarArrow(m.Gram(a.Inverse().Transpose())) << "\n"
        ExportVarArrow(m) << "\n"
        ExportVarArrow(a) << "\n"
        << "\n";
        //        m0 = m;
        //       std::cout ExportVarArrow
    }
    */
    
    // Sort and reorder
    for(int i=d; i>=1; --i){
        const int j=i-1;
        // Break if the norms are correctly sorted
        if(m(i,i) >= m(j,j)) return i+1;
        
        // Otherwise exchange the two vectors
        VectorType & rowi = *(VectorType*)&a(i,0);
        VectorType & rowj = *(VectorType*)&a(j,0);
        std::swap(rowi, rowj);

        // Exchange the scalar products as well
        std::swap(m(i,i), m(j,j));
        for(int k=0;   k<j; ++k) {std::swap(m(i,k), m(j,k));}
        for(int k=i+1; k<Dimension; ++k) {std::swap(m(i,k), m(j,k));}
    }
    return 1;
}

template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
IsMinkowskiReduced(const SymmetricMatrixType & m) -> ScalarType {
    ScalarType tol = 0;
    for(int d=1; d<Dimension; ++d){
        const ScalarType normD2 = m(d,d);
        tol = std::max(tol, m(d-1,d-1)/normD2 ); // Check that sorted
        VectorType c = VectorType::Constant(0);
        for(c[d]=1; c[d]< (d==5 ? 3 : 2); ++c[d]){
            for(int k=1; k< (1<< (2*d)); ++k){
                bool skip=false;
                for(int j=0; j<d; ++j){
                    switch((k >> (2*j)) & 3){
                        case 0: c[j]=0; break;
                        case 1: c[j]=1; break;
                        case 2: c[j]=-1; break;
                        case 3: skip=true;
                    }
                }
                if(skip) continue; // Possible improvement: increase k suitably ?
/*                if(normD2/m.SquaredNorm(c)>1){
                    std::cout
                    ExportVarArrow(c)
                    ExportVarArrow(m)
                    << "\n";
                }*/
                tol = std::max(tol, normD2/m.SquaredNorm(c) ); // Check that cannot be reduced
            } // for k
        } // for c[d]
    } // for d
    return tol-1;
}

#endif /* VoronoiReduction_Basis_hxx */
