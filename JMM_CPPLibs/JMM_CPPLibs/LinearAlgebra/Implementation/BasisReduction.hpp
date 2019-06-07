// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef BasisReduction_hpp
#define BasisReduction_hpp

template<typename TS, typename TD, size_t VD> const typename BasisReduction<TS,TD,VD>::BasisType
BasisReduction<TS,TD,VD>::CanonicalBasis() {BasisType basis;
    for(int i=0; i<Dimension; ++i) for(int j=0; j<Dimension; ++j) basis[i][j]=(i==j);
    return basis;}

template<typename TS, typename TD, size_t VD> const typename BasisReduction<TS,TD,VD>::SuperbaseType
BasisReduction<TS,TD,VD>::CanonicalSuperBase() {SuperbaseType sb;
    for(int i=0; i<=Dimension; ++i) for(int j=0; j<Dimension; ++j) sb[i][j]= i<Dimension ? (i==j) : -1;
    return sb;}

template<typename TS, typename TD, size_t VD> void
BasisReduction<TS,TD,VD>::ReducedBasis(const SymmetricMatrixType & m, BasisType & basis) {
//    std::cout << "\nReducing " << m << "\n";
    if(Dimension==1) return;
    std::array<ScalarType,Dimension> sqNorms;
    for(int i=0; i<Dimension; ++i) sqNorms[i]=m.SquaredNorm(basis[i]);
    
    std::array<int,Dimension> order;
    for(int i=0; i<Dimension; ++i) order[i]=i;
    std::sort(order.begin(), order.end(), [&](int i,int j){return sqNorms[i]<sqNorms[j];});
    
    auto b = [&](int i)->DiscreteVectorType & {return basis[order[i]];};
    auto sqn = [&](int i)-> ScalarType & {return sqNorms[order[i]];};

//    std::cout ExportArrayArrow(basis)  ExportArrayArrow(sqNorms) ExportArrayArrow(order) << "\n";

    while(true){
        const ScalarType k01 = m.ScalarProduct(b(0), b(1)) / sqn(0);
        if(std::abs(k01) > 1./2){
            b(1)-=std::round(k01)*b(0);
            sqn(1)=m.SquaredNorm(b(1));
            if(Dimension==2 || sqn(0)>sqn(1)) std::swap(order[0],order[1]);
//            std::cout ExportArrayArrow(basis)  ExportArrayArrow(sqNorms) ExportArrayArrow(order) << "\n";
            continue;
        } else if(Dimension==2) break;

        if(Dimension==3){
            const ScalarType
            k02 = m.ScalarProduct(b(0),b(2))/sqn(0),
            k12 = m.ScalarProduct(b(1),b(2))/sqn(1);
            
            if(std::abs(k02)>1./2.){
                b(2)-=std::round(k02)*b(0);
                sqn(2)=m.SquaredNorm(b(2));
            } else if(std::abs(k12)>1./2.){
                b(2)-=std::round(k12)*b(1);
                sqn(2)=m.SquaredNorm(b(2));
            } else break;
            
            if(sqn(1)>sqn(2)) {
                std::swap(order[1],order[2]);
                if(sqn(0)>sqn(1)) std::swap(order[0],order[1]);
            }
            continue;
        } // Dimension==3
    } // while
}

template<typename TS, typename TD, size_t VD> void
BasisReduction<TS,TD,VD>::ObtuseSuperbase(const SymmetricMatrixType & m, SuperbaseType & sb) {
    static_assert(VD<=3,"Sorry, dimensions >=4 are not supported");
    if(m.Determinant()<=0 || m.Trace()<=0)
        ExceptionMacro("ObtuseSuperbase error : input must be positive definite");
    if(Dimension==1) return;
    bool reduced=true;
    for(int countIt=0; reduced && countIt<maxIt; ++countIt){
        reduced=false;
        for(int i=0; i<Dimension; ++i)
            for(int j=i+1; j<=Dimension; ++j){
                if(m.ScalarProduct(sb[i],sb[j])>0){
                    reduced=true;
                    if(Dimension==2){
                        const int k = 0+1+2-i-j;
                        sb[k]=sb[i]-sb[j];
                    } else {
                        int k[2]; int *kIt=k;
                        for(int l=0; l<=Dimension; ++l)
                            if(l!=i && l!=j){
                                *kIt=l;
                                ++kIt;
                            }
                        assert(kIt==k+2);
                        assert(i+j+k[0]+k[1]==0+1+2+3);
                        sb[k[0]]+=sb[i];
                        sb[k[1]]+=sb[i];
                    }
                    sb[i]*=-1;
                }
            }
    }
    if(reduced) assert(false);
    if(reduced) ExceptionMacro("ObtuseSuperbase error : did not terminate in "
                               << maxIt << " iterations, for matrix " << m
                               /*<< ", current superbase: " ExportArrayArrow(sb)*/
                               << ".\n");
    for(int i=0; i<Dimension; ++i)
        for(int j=i+1; j<=Dimension; ++j)
            assert(m.ScalarProduct(sb[i],sb[j])<=0);
}


template<typename TS, typename TD, size_t VD> typename BasisReduction<TS,TD,VD>::TensorDecompositionType
BasisReduction<TS,TD,VD>::TensorDecomposition(const SymmetricMatrixType & diff){
    SuperbaseType sb=CanonicalSuperBase();
    ObtuseSuperbase(diff, sb);
    return TensorDecompositionHelper::Get(diff,sb);
}

template<typename TS, typename TD, size_t VD> template<typename Dummy>
struct BasisReduction<TS, TD, VD>::TensorDecompositionHelper_<1,Dummy> {
    static_assert(VD==1,"Inconsistent dimensions");
    typedef BasisReduction<TS, TD, VD> T;
	Redeclare3Types(T, TensorDecompositionType, SymmetricMatrixType, SuperbaseType);
    static TensorDecompositionType Get(const SymmetricMatrixType & diff, const SuperbaseType &){
        TensorDecompositionType decomp;
        decomp.offsets[0][0]=1;
        decomp.weights[0] = diff(0,0);
        return decomp;
    }
};
    
template<typename TS, typename TD, size_t VD> template<typename Dummy>
struct BasisReduction<TS, TD, VD>::TensorDecompositionHelper_<2,Dummy> {
    static_assert(VD==2,"Inconsistent dimensions");
    typedef BasisReduction<TS, TD, VD> T;
	Redeclare3Types(T, TensorDecompositionType, SymmetricMatrixType, SuperbaseType);
    static TensorDecompositionType Get(const SymmetricMatrixType & diff, const SuperbaseType & sb){
        TensorDecompositionType decomp;
        for(int i=0; i<3; ++i){
            int j=(i+1)%3,k=(i+2)%3;
            decomp.offsets[i] = LinearAlgebra::Perp(sb[i]);
            decomp.weights[i] = -diff.ScalarProduct(sb[j],sb[k]);
        }
        return decomp;
    }
};

template<typename TS, typename TD, size_t VD> template<typename Dummy>
struct BasisReduction<TS, TD, VD>::TensorDecompositionHelper_<3,Dummy> {
    static_assert(VD==3,"Inconsistent dimensions");
    typedef BasisReduction<TS, TD, VD> T;
	Redeclare3Types(T, TensorDecompositionType, SymmetricMatrixType, SuperbaseType);
    static TensorDecompositionType Get(const SymmetricMatrixType & diff, const SuperbaseType & sb){
        TensorDecompositionType decomp;
        auto offsetIt = decomp.offsets.begin();
        auto weightIt = decomp.weights.begin();
        for(int i=0; i<4; ++i)
            for(int j=i+1; j<4; ++j){
                int rem[2]; int*remIt=rem;
                for(int k=0; k<4; ++k) if(k!=i && k!=j) {*remIt=k; ++remIt;}
                *offsetIt = LinearAlgebra::Cross(sb[rem[0]],sb[rem[1]]);
                *weightIt = -diff.ScalarProduct(sb[i],sb[j]);
                ++offsetIt; ++weightIt;
            }
        return decomp;
    }
};


#endif /* BasisReduction_hpp */
