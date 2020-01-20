// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_LinearTransform_h
#define AmongConvex2_LinearTransform_h

#include "VectorType.h"
#include "../DataStructures/GetComponent.h"

namespace LinearAlgebra {
    
template<typename TComponent, size_t VRows, size_t VColumns>
struct Matrix :
vector_space< Matrix<TComponent,VRows,VColumns>, TComponent>
{
    typedef TComponent ComponentType;
    static constexpr size_t Rows = VRows, Columns = VColumns;
    
    static constexpr bool IsInRange(int i, int j)  {
		return 0<=i && i<Rows && 0<=j && j<Columns;}
    ComponentType           & operator()(size_t i, size_t j)       {
        assert(IsInRange(int(i),int(j))); return data[LinearizedIndex(i,j)];}
    const ComponentType     & operator()(size_t i, size_t j) const {
        assert(IsInRange(int(i),int(j))); return data[LinearizedIndex(i,j)];}
    
    // Linear operations
    Matrix & operator+=(const Matrix & m){data+=m.data; return *this;}
    Matrix & operator-=(const Matrix & m){data-=m.data; return *this;}
    Matrix & operator*=(const ComponentType & a){data*=a;   return *this;}
    Matrix & operator/=(const ComponentType & a){data/=a;   return *this;}
    Matrix operator -() {Matrix m; m.data = -data; return m;}
    
    // Multiplication
    typedef Vector<ComponentType,Rows>    OutputVectorType;
    typedef Vector<ComponentType,Columns>  InputVectorType;
	template<typename T, typename S = algebra_t<ComponentType, T> >
    Vector<S,Rows> operator * (const Vector<T, Columns> & u) const;
	OutputVectorType Column(int) const;
	InputVectorType  Row(int) const;
    
    template <size_t Columns2> Matrix<ComponentType, Rows, Columns2>
    operator * (const Matrix<ComponentType,Columns,Columns2> & m) const;
    
    // Matrix algebra
    ComponentType Trace() const;
    ComponentType Determinant() const;
    Matrix<ComponentType,Columns,Rows> Transpose() const;

    ComponentType FrobeniusSquaredNorm() const {return data.SquaredNorm();}
    ComponentType FrobeniusNorm() const {return data.Norm();}
    
    Matrix Inverse() const;
    InputVectorType Solve(OutputVectorType) const;
    
    // Constructors
    void fill(const ComponentType & a){data.fill(a);}
	static Matrix Zero(){Matrix a; a.fill(ComponentType(0)); return a;}
    static Matrix Identity();
    // Do not forget double braces, like {{ {{1,2}} , {{3,4}} }}.
    static Matrix FromRows(const std::array<InputVectorType,Rows> & t);
    static Matrix FromColumns(const std::array<OutputVectorType,Columns> & t);
    static Matrix Rotation(ComponentType theta); // 2x2 only
    static Matrix Random(); // random entries uniform in [-1,1]
    template<typename T> static Matrix CastCoordinates(const Matrix<T,Rows,Columns> &);

    static const size_t InternalDimension = Rows*Columns;
	Vector<ComponentType, InternalDimension> data;
	template<typename ...T,typename dummy = typename std::enable_if<sizeof...(T)==InternalDimension>::type >
    constexpr Matrix(T... t):data(t...){};

    Matrix(){};
protected:
    static size_t LinearizedIndex(size_t i, size_t j) {return i*Columns+j;} // Row major
//    template<size_t d, typename Dummy> struct Det;
};

// Row based printing.
template<typename TC, size_t VR, size_t VC>
std::ostream & operator << (std::ostream & f, const Matrix<TC,VR,VC> & m) {
    f<<"{";
    for(int i=0; i<VR; ++i){
        if(i>0) f<<",";
        f<<"{";
        for(int j=0; j<VC; ++j){
            if(j>0) f<<",";
            f<<m(i,j);
        }
        f<<"}";
    }
    f<<"}";
    return f;
}
    
    
#include "Implementation/MatrixType.hxx"

} // namespace LinearAlgebra

template<typename C, size_t VRows, size_t VCols>
struct GetComponent<LinearAlgebra::Matrix<C,VRows,VCols>, C> {
	typedef LinearAlgebra::Matrix<C,VRows,VCols> T;
	static constexpr size_t size() {return T().data.size();}
	static const C & Get(const T & t, size_t i) {assert(i<size()); return t.data[i];}
	static C & Get(T & t, size_t i) {assert(i<size()); return t.data[i];}
};

#endif
