//
//  HopfLaxMinimize.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 21/08/2018.
//

#ifndef HopfLaxMinimize_h
#define HopfLaxMinimize_h

/**
 This file implements methods for minimizing
 F(u) + \<d,u\> on the set {\<1,u\>=1, u>0}
 where F is a norm which is of either Riemann, Rander, or asymmetric quadratic type.
 The minimiser u is returned as well. (Akin to a covariant gradient.)
 
 Caution : the value +Infinity is returned if the infimum is not attained on the
 relative interior of the set {\<1,u\>=1, u>0}
 
 */

// ****** Minimize . Infinity if not attained in interior.*****

#include "SymmetricMatrixType.h"
#include "SymmetricMatrixPair.h"
#include "RanderNorm.h"
#include "AsymmetricQuadraticNorm.h"

namespace LinearAlgebra {
    
    template<typename TC, size_t VD>
    TC CoefficientsSum(const SymmetricMatrix<TC, VD> & m){
        TC result(0);
        for(int i=0; i<VD; ++i)
            for(int j=0; j<=i; ++j)
                result+= (i==j ? m(i,j) : 2*m(i,j));
        return result;
    }
    
    template<typename TS0, typename TS1>
    typename TS0::ComponentType CoefficientsSum(const SymmetricMatrixPair<TS0, TS1> & m){
        return CoefficientsSum(m.m0)+CoefficientsSum(m.m1);
    }
    
    template<
    typename SymmetricMatrixType,
    typename VectorType = typename SymmetricMatrixType::VectorType,
    typename ComponentType = typename SymmetricMatrixType::ComponentType
    >
    ComponentType
    HopfLaxMinimizeInv_NoCheck(const SymmetricMatrixType & m,
                               const VectorType & d,
                               VectorType & g){
		using std::sqrt;
		if(VectorType::Dimension==1){
			const ComponentType invNorm0 = sqrt(m(0,0));
			g[0]=invNorm0;
			return 1/invNorm0 + d[0];
		}
		
        VectorType ones;
        ones.fill(1);
        
        // Introduce the quadratic system a x^2 - 2 b x + c, for |lambda*ones-d|^2=1
        /* // We use slightly optimized equivalents of the following.
         const ComponentType
         a = m.SquaredNorm(ones),
         b = m.ScalarProduct(ones, d),
         c = m.SquaredNorm(d) - 1;
         */
        const ComponentType a = CoefficientsSum(m);
        const VectorType md = m*d;
        const ComponentType
        b = md.Sum(),
        c = d.ScalarProduct(md) -1;
        
        const ComponentType delta2 = b*b-a*c;
        if(delta2<0) return std::numeric_limits<ComponentType>::infinity();
        
        const ComponentType delta = sqrt(delta2);
        const ComponentType result = (b + delta)/a; // Largest root
        
        g = m*(result*ones-d);
        if( !g.IsPositive() ) return std::numeric_limits<ComponentType>::infinity();
        
        // g/=g.Sum(); // Now return the geodesic flow, not barycentric coordinates
        return result;
    }
    
    
    template<
    typename SymmetricMatrixType,
    typename VectorType = typename SymmetricMatrixType::VectorType,
    typename ComponentType = typename SymmetricMatrixType::ComponentType
    >
    ComponentType
    HopfLaxMinimize(const SymmetricMatrixType & gram,
                    const VectorType & d,
                    VectorType & g){
        
        if(!d.IsFinite()) return std::numeric_limits<ComponentType>::infinity();
        else return HopfLaxMinimizeInv_NoCheck(gram.Inverse(), d,g);
    }
    
    
    template<
    typename SymmetricMatrixType,
    typename VectorType = typename SymmetricMatrixType::VectorType,
    typename ComponentType = typename SymmetricMatrixType::ComponentType
    >
    ComponentType
    HopfLaxMinimizeInv(const SymmetricMatrixType & gramInv,
                       const VectorType & d,
                       VectorType & g){
        
        if(!d.IsFinite()) return std::numeric_limits<ComponentType>::infinity();
        else return HopfLaxMinimizeInv_NoCheck(gramInv, d,g);
    }
    
    template<typename ScalarType, size_t Dimension>
    ScalarType HopfLaxMinimize(const RanderNorm<ScalarType,Dimension> & norm,
                               const Vector<ScalarType,Dimension> & d,
                               Vector<ScalarType,Dimension> & g){
		const ScalarType val = HopfLaxMinimize(norm.m, d-norm.w, g);
		g/=(1-norm.w.ScalarProduct(g)); // Normalize the gradient
		return val;
    }
    
    template<typename ScalarType, size_t Dimension>
    ScalarType HopfLaxMinimize(const AsymmetricQuadraticNorm<ScalarType, Dimension> & norm,
                               const Vector<ScalarType, Dimension> & d,
                               Vector<ScalarType, Dimension> & g){
        if(Dimension==1){
			const ScalarType norm0 = sqrt(norm.m(0,0)+ square(std::min(0.,norm.w[0])));
            g[0]=1/norm0;
            return norm0 + d[0];
        }
		
        if( 0 <= *std::max_element(norm.w.begin(), norm.w.end()) ){
            const ScalarType val = HopfLaxMinimize(norm.m, d, g);
            const ScalarType wMean = g.ScalarProduct(norm.w);
            if(wMean>=0) return val;
        }
        using SymmetricMatrixType = SymmetricMatrix<ScalarType, Dimension>;
        const ScalarType val = HopfLaxMinimize(norm.m+SymmetricMatrixType::RankOneTensor(norm.w), d, g);
		
		if(g.ScalarProduct(norm.w)>0) {return std::numeric_limits<ScalarType>::infinity();}
        return val;
        }
    }

#endif /* HopfLaxMinimize_h */
