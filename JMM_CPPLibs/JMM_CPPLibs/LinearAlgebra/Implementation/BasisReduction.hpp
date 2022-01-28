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
	const bool invalid = m.Determinant()<=0 || m.Trace()<=0;
	assert(!invalid); // Detects (some, not all) non positive definite matrices.
    if(invalid) ExceptionMacro("ObtuseSuperbase error : input must be positive definite");
    if(Dimension==1) return;
    bool reduced=true;
    for(int countIt=0; reduced && countIt<maxIt; ++countIt){
        reduced=false;
        for(int i=0; i<Dimension; ++i)
            for(int j=i+1; j<=Dimension; ++j){
                if(m.ScalarProduct(sb[i],sb[j])>0){
                    reduced=true;
                    if constexpr(Dimension==2){
						const auto [k] = ComplementIndices(i,j);
                        sb[k]=sb[i]-sb[j];
                    } else {
						const auto [k,l] = ComplementIndices(i,j);
                        sb[k]+=sb[i];
                        sb[l]+=sb[i];
                    }
                    sb[i]*=-1;
                }
            }
    }
    if(reduced) assert(false);
    if(reduced) ExceptionMacro("ObtuseSuperbase error : did not terminate in "
                               << maxIt << " iterations, for matrix " << m
                               // << ", current superbase: " ExportArrayArrow(sb)
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


template<typename TS, typename TD, size_t VD> auto BasisReduction<TS, TD, VD>::
ComplementIndices(int i, int j) -> std::array<int,Dimension-1> {
	assert(0<=i && i<=Dimension);
	assert(0<=j && j<=Dimension);

	if constexpr(Dimension==2){return {(0+1+2)-i-j};}
	else { static_assert(Dimension==3, "Unsupported dimension");
		std::array<int,2> result;
		auto rIt = result.begin();
		for(int k=0; k<=Dimension; ++k){
			if(k!=i && k!=j){
				*rIt = k;
				++rIt;
			}
		}
		assert(rIt==result.end());
		assert(i+j+result[0]+result[1]==0+1+2+3);
		return result;
	}
}

template<typename TS, typename TD, size_t VD> int BasisReduction<TS, TD, VD>::
LinearizeIndices(int i, int j){
	if(i>j) std::swap(i,j);
	assert(0<=i && i<j && j<=Dimension);
	int r=j-i-1;
	if(i==1){r+=Dimension;}
	else if(i==2 && Dimension==3){r+=2*Dimension-1;}
	else {assert(i==0);}
	assert(r == i*Dimension-(i*(i-1))/2 + (j-i-1) );
	return r;
}

// ------------ Selling path -----------

template<typename TS,typename TD, size_t VD> void BasisReduction<TS,TD,VD>::SellingPath::
PrintSelf(std::ostream & os) const {
	os << "{"
	ExportVarArrow(D0)
	ExportVarArrow(D1)
	ExportArrayArrow(sb)
	ExportArrayArrow(offsets)
	ExportArrayArrow(weights0)
	ExportArrayArrow(weights1)
	<< "}";
}

template<typename TS, typename TD, size_t VD> BasisReduction<TS, TD, VD>::SellingPath::
SellingPath(const SymmetricMatrixType & D0_, const SymmetricMatrixType & D1_, ScalarType t)
:D0(D0_),D1(D1_),sb(CanonicalSuperBase()){
	
	// Find the obtuse superbase for t.
	ObtuseSuperbase((1.-t)*D0 + t*D1, sb);
	
	// Fill in the weights and offsets
	for(int i=0, r=0; i<=Dimension; ++i){
		for(int j=i+1; j<=Dimension; ++j){
			weights0[r] = -D0.ScalarProduct(sb[i],sb[j]);
			weights1[r] = -D1.ScalarProduct(sb[i],sb[j]);
			if constexpr(Dimension==2){
				const auto [k] = ComplementIndices(i,j);
				offsets[r] = Perp(sb[k]);
			} else {static_assert(Dimension==3,"Unsupported dimension");
				const auto [k,l] = ComplementIndices(i,j);
				offsets[r] = Cross(sb[k],sb[l]);
			}
			++r;
		}
	}
}

template<typename TS, typename TD, size_t VD> auto BasisReduction<TS,TD,VD>::SellingPath::
NextStep(ScalarType & t) const -> std::pair<int,int> {
	[[maybe_unused]] const ScalarType tMin = t;
	[[maybe_unused]] constexpr ScalarType eps = std::numeric_limits<ScalarType>::epsilon();

	t=std::numeric_limits<ScalarType>::infinity();
	// Find the next t, when obtuseness fails.
	int iNext=-1,jNext=-1;
	for(int i=0,r=0; i<=Dimension; ++i) {
		for(int j=i+1; j<=Dimension; ++j,++r) {
			const ScalarType s0 = weights0[r], s1=weights1[r];
//			std::cout ExportVarArrow((1.-tMin)*s0+tMin*s1) << std::endl;
			assert((1.-tMin)*s0+tMin*s1 >= -1e-11); //Expect obtuse superbase initially
			if(s1>=s0) continue; // Continue if vectors become more obtuse
			if(0>=s0) continue;  // Only useful in degenerate cases
			const ScalarType tRoot = s0/(s0-s1);
			assert(tRoot>=tMin-1000*eps); //Expect obtuse superbase initially
			if(tRoot>=t) continue;
			t=tRoot;
			iNext=i;
			jNext=j;
		}
	}
	return {iNext,jNext};
}

template<typename TS, typename TD, size_t VD> void BasisReduction<TS,TD,VD>::SellingPath::
MakeStep(const std::pair<int,int> & indices) {
	// Make a Selling step
	const auto [i,j] = indices;
	
	if constexpr(Dimension==2){
		const auto [k] = ComplementIndices(i,j);
		// Update superbase (Selling step)
		sb[k]=sb[i]-sb[j];
		sb[i]*=-1;
		// Update offsets (sign changes omitted)
		offsets[LinearizeIndices(i,j)] = Perp(sb[k]);
	} else { static_assert(Dimension==3, "Unsupported dimension");
		const auto [k,l] = ComplementIndices(i,j);
		// Update superbase (Selling step)
		sb[k]+=sb[i];
		sb[l]+=sb[i];
		sb[i]*=-1;
		// Update offsets (sign changes omitted)
		offsets[LinearizeIndices(i,j)] = Cross(sb[k],sb[l]);
		std::swap(offsets[LinearizeIndices(i, k)],offsets[LinearizeIndices(i, l)]);
	}
	
	// Update weights
	for(int i=0,r=0; i<=Dimension; ++i){
		for(int j=i+1; j<=Dimension; ++j,++r){
			weights0[r] = -D0.ScalarProduct(sb[i],sb[j]);
			weights1[r] = -D1.ScalarProduct(sb[i],sb[j]);
		}
	}
}

#endif /* BasisReduction_hpp */
