// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef LiftedFastMarching_SymmetricMatrixType_h
#define LiftedFastMarching_SymmetricMatrixType_h

#include "MatrixType.h"
#include "../DataStructures/GetComponent.h"
#include "../Macros/Exception.h"

namespace LinearAlgebra {
    
template<typename TComponent, size_t VDimension>
struct SymmetricMatrix :
vector_space< SymmetricMatrix<TComponent, VDimension>, TComponent>
{
    typedef TComponent ComponentType;
    static constexpr size_t Dimension = VDimension;
    typedef Vector<ComponentType, Dimension> VectorType;
    typedef Point<ComponentType, Dimension> PointType;
    typedef Matrix<ComponentType, Dimension, Dimension> MatrixType;
    
    // Coefficient access
    static constexpr bool IsInRange(int i, int j){
        return 0<=i && i<Dimension && 0<=j && j<Dimension;}
    static constexpr int LinearizedIndex(int i, int j);
    ComponentType & operator()(int i, int j){return data[LinearizedIndex(i,j)];}
    const ComponentType & operator()(int i, int j) const {return data[LinearizedIndex(i,j)];}

    // Linear operations
    SymmetricMatrix & operator+=(const SymmetricMatrix & m){data+=m.data; return *this;}
    SymmetricMatrix & operator-=(const SymmetricMatrix & m){data-=m.data; return *this;}
    SymmetricMatrix & operator*=(const ComponentType & a){data*=a; return *this;}
    SymmetricMatrix & operator/=(const ComponentType & a){data/=a; return *this;}
    SymmetricMatrix operator -() const {SymmetricMatrix m; m.data=-data; return m;}
    
    // Geometry
    template<typename T,typename S = algebra_t<ComponentType, T> > S
    ScalarProduct(const Vector<T,Dimension> &, const Vector<T, Dimension> &) const;
    
    template<typename T,typename S = algebra_t<ComponentType, T> > S
    SquaredNorm(const Vector<T,Dimension> & u) const {return ScalarProduct(u, u);}
    
    template<typename T, typename S = algebra_t<ComponentType, T> > S
	Norm(const Vector<T,Dimension> & u) const {using std::sqrt;return sqrt(SquaredNorm(u));}
    
    template<typename T> bool
    IsAcute(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v) const {
        return ScalarProduct(u, v)>=0;}

    template<typename T, typename S = algebra_t<T,ComponentType> > Vector<S,Dimension>
    operator*(const Vector<T, Dimension> &) const;
	
	VectorType Gradient(const VectorType & u) const {return operator*(u)/Norm(u);}
	
    ComponentType Trace() const;
    ComponentType Determinant() const;
	SymmetricMatrix Comatrix() const;
    SymmetricMatrix Inverse() const;
    VectorType CGSolve(VectorType) const;
	
	// Cholevski factorisation of positive matrices
	typedef std::pair<bool,MatrixType> CholevskiType;
	CholevskiType Cholevski() const;
	bool IsPositive() const {return Cholevski().first;}
	
	
/*	// Spectral norm by the power method
	SymmetricMatrix Square() const;
	ComponentType SpectralNorm() const;*/
	
	ComponentType SquaredFrobeniusNorm() const;
	ComponentType FrobeniusNorm() const;
    
    // Factory methods
    static SymmetricMatrix Zero();
    static SymmetricMatrix Identity();
    static SymmetricMatrix Diagonal(const VectorType &);
	/// Returns v v^T
	template<typename T> static SymmetricMatrix OuterSelf(const Vector<T,Dimension> & u);
    template<typename T> static SymmetricMatrix
	RankOneTensor(const Vector<T, Dimension> & u) {return OuterSelf(u);}
	/// Returns u v^T+v u^T
	static SymmetricMatrix Outer2(const VectorType &u, const VectorType &v);
	
    static SymmetricMatrix RandomPositive();
    static SymmetricMatrix FromUpperTriangle(const MatrixType & mat);
    SymmetricMatrix operator=(ComponentType a); // Dimension 1 only
    
	/// Returns symmetric matrix M_ij = m.ScalarProduct(a[i],a[j]);*/
    template<typename T, size_t D> SymmetricMatrix<ComponentType,D>
    Gram(const std::array<Vector<T,Dimension>, D> & a) const;
    template<size_t D> SymmetricMatrix<ComponentType, D> /// a^T.m.a
    Gram(const Matrix<ComponentType,Dimension,D> & a) const;
    template<size_t D> SymmetricMatrix<ComponentType, D> /// a.m.a^T
    GramT(const Matrix<ComponentType,D,Dimension> & a) const;
    template<typename T, size_t D> static SymmetricMatrix
	EuclideanGram(const std::array<Vector<T,D>, Dimension> & a);

    // Conversion
    explicit operator MatrixType() const;
    template<typename T> static SymmetricMatrix CastCoordinates(const SymmetricMatrix<T, Dimension> & m0);
    
    static const size_t InternalDimension = (Dimension*(Dimension+1))/2;
    Vector<ComponentType,InternalDimension> data;
    template<typename ...T,typename dummy = typename std::enable_if<sizeof...(T)==InternalDimension>::type >
    constexpr SymmetricMatrix(T... t):data(t...){};
	
	SymmetricMatrix ComponentWiseProduct(const SymmetricMatrix & other) const {
		SymmetricMatrix res; res.data=data.ComponentWiseProduct(other.data); return res;}
	
    SymmetricMatrix(){};
protected:
    const ComponentType & coef(int i, int j) const {return this->operator()(i,j);}
};


    // Printing
    template<typename TC, size_t VD>
    std::ostream & operator << (std::ostream & f, const SymmetricMatrix<TC,VD> & m)
    {
        f<<"{";
        for(int i=0; i<VD; ++i){
            if(i>0) f<<",";
            f<<"{";
            for(int j=0; j<VD; ++j){
                if(j>0) f<<",";
                f<<m(i,j);
            }
            f<<"}";
        }
        f<<"}";
        return f;
    }

#include "Implementation/SymmetricMatrixType.hxx"
} // end namespace LinearAlgebra

template<typename C, size_t VD>
struct GetComponent<LinearAlgebra::SymmetricMatrix<C,VD>, C> {
	typedef LinearAlgebra::SymmetricMatrix<C,VD> T;
	static constexpr size_t size() {return T().data.size();}
	static const C & Get(const T & t, size_t i) {assert(i<size()); return t.data[i];}
	static C & Get(T & t, size_t i) {assert(i<size()); return t.data[i];}
};

#endif
