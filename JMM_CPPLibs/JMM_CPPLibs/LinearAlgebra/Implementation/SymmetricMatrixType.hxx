//
//  SymmetricMatrixType.hxx
//  
//
//  Created by Jean-Marie Mirebeau on 09/02/2018.
//

#ifndef SymmetricMatrixType_hxx
#define SymmetricMatrixType_hxx

template<typename TC, size_t VD> constexpr int SymmetricMatrix<TC,VD>::
LinearizedIndex(int i, int j){
    assert(IsInRange(i, j));
    if(i<j) std::swap(i, j);
    return (i*(i+1))/2+j;
}

template<typename TC, size_t VD> template<typename T, typename S>
S SymmetricMatrix<TC,VD>::
ScalarProduct(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v) const
{
    S sum(0.);
    for(int i=0; i<Dimension; ++i){
        S sumi(0.);
		for(int j=0; j<Dimension; ++j){
			sumi+=coef(i,j)*v[j];}
		sum+=u[i]*sumi;
	}
    return sum;
}

template<typename TC, size_t VD> template<typename T, typename S>
auto SymmetricMatrix<TC,VD>::
operator*(const Vector<T,Dimension> & u) const -> Vector<S,Dimension>
{
    Vector<S,Dimension> v;
    for(int i=0; i<Dimension; ++i){
        v[i]=S(0.);
		for(int j=0; j<Dimension; ++j){
			v[i]+=coef(i,j)*u[j];
		}
	}
    return v;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Zero() -> SymmetricMatrix
{
    SymmetricMatrix m;
    std::fill(m.data.begin(), m.data.end(), ComponentType(0));
    return m;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Identity() -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)=(i==j);
    return m;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Diagonal(const VectorType & u) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)= i==j ? u[i] : ComponentType(0);
    return m;
}

template<typename TC, size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
OuterSelf(const Vector<T,Dimension> & u) -> SymmetricMatrix {
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)=u[i]*u[j];
    return m;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Outer2(const VectorType & u, const VectorType & v) -> SymmetricMatrix {
	SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)=u[i]*v[j] + u[j]*v[i];
    return m;
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
RandomPositive() -> SymmetricMatrix {    
    MatrixType a=MatrixType::Random();
    return FromUpperTriangle(a.Transpose()*a);
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
FromUpperTriangle(const MatrixType & mat) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=i; j<Dimension; ++j)
            m(i,j)=mat(i,j);
            return m;
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
operator=(ComponentType a) -> SymmetricMatrix
{
    static_assert(Dimension==1,"Assignment from scalar to matrix in dimension 1 only");
    data[0]=a;
    return *this;
}


template<typename TC, size_t VD> template<typename T, size_t D> auto SymmetricMatrix<TC,VD>::
Gram(const std::array<Vector<T,Dimension>, D> & a) const
-> SymmetricMatrix<ComponentType,D> {
    SymmetricMatrix<ComponentType,D> m;
    for(int i=0; i<D; ++i){
        const Vector<ComponentType, Dimension> mai = operator*(a[i]);
        for(int j=0; j<=i; ++j)
            m(i,j) = mai.ScalarProduct(a[j]);
    }
    /*
     SymmetricMatrix<ComponentType,D> m;
     for(int i=0; i<D; ++i)
     for(int j=0; j<=i; ++j)
     m(i,j) = ScalarProduct(a[i],a[j]);*/
    return m;
}

template<typename TC, size_t VD> template<size_t D> auto SymmetricMatrix<TC,VD>::
Gram(const Matrix<ComponentType,Dimension,D> & a) const
-> SymmetricMatrix<ComponentType, D>
{
	std::array<Vector<ComponentType,Dimension>,D> arr;
	for(int i=0; i<D; ++i){
		for(int j=0; j<Dimension; ++j){
			arr[i][j] = a(j,i);
		}
	}
	return Gram(arr);
}

template<typename TC, size_t VD> template<size_t D> auto SymmetricMatrix<TC,VD>::
GramT(const Matrix<ComponentType,D,Dimension> & a) const
-> SymmetricMatrix<ComponentType, D>
{
	std::array<Vector<ComponentType,Dimension>,D> arr;
	for(int i=0; i<D; ++i){
		for(int j=0; j<Dimension; ++j){
			arr[i][j] = a(i,j);
		}
	}
	return Gram(arr);
}

template<typename TC, size_t VD> template<typename T, size_t D> auto SymmetricMatrix<TC,VD>::
EuclideanGram(const std::array<Vector<T,D>, Dimension> & a) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j) = a[i].ScalarProduct(a[j]);
    return m;
}

// Conversion
template<typename TC,size_t VD> SymmetricMatrix<TC,VD>::
operator MatrixType() const {
    MatrixType result;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<Dimension; ++j)
            result(i,j)=coef(i, j);
    return result;
}

template<typename TC,size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
CastCoordinates(const SymmetricMatrix<T, Dimension> & m0) -> SymmetricMatrix {
    SymmetricMatrix m;
    m.data = Vector<ComponentType,InternalDimension>::CastCoordinates(m0.data);
    return m;
}


// Matrix operations
template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
Trace() const -> ComponentType
{
    ComponentType result=0;
    for(int i=0; i<Dimension; ++i)
        result+=coef(i,i);
        return result;
}


template<typename TC, size_t VD>
TC
SymmetricMatrix<TC, VD>::Determinant() const {
	if constexpr(VD==0) {return ComponentType(1.);}
	else if constexpr(VD==1) {return coef(0,0);}
	else if constexpr(VD==2) {return coef(0,0)*coef(1,1)-coef(1,0)*coef(0,1);}
	else if constexpr(VD==3) {return
		coef(0,0)*(coef(1,1)*coef(2,2)-coef(1,2)*coef(2,1))
		+coef(0,1)*(ComponentType(2.)*coef(1,2)*coef(2,0)-coef(1,0)*coef(2,2))
		-coef(0,2)*(coef(1,1)*coef(2,0));
}
	else {
		LinearAlgebra::Matrix<TC, VD, VD> mat;
		for(int i=0; i<VD; ++i)
			for(int j=0; j<VD; ++j)
				mat(i,j)=coef(i,j);
		return mat.Determinant();
	}
	
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC, VD>::
Comatrix() const -> SymmetricMatrix {
	static_assert(Dimension<=3, "Unsupported dimension");
	if constexpr(Dimension==0){return SymmetricMatrix();}
	else if constexpr(Dimension==1){return SymmetricMatrix{1.};}
	else if constexpr(Dimension==2){
		return SymmetricMatrix{coef(1,1),-coef(1,0),coef(0,0)};}
	else {
		static_assert(Dimension==3, "Unsupported dimension");
		SymmetricMatrix  m;
		for(int i=0; i<3; ++i)
			for(int j=0; j<=i; ++j){
				const int i1 = (i+1)%3, i2=(i+2)%3, j1=(j+1)%3, j2=(j+2)%3;
				m(i,j) = coef(i1,j1)*coef(i2,j2)-coef(i1,j2)*coef(i2,j1);
			}
		return m;
	}
}

template<typename TC, size_t VD> SymmetricMatrix<TC, VD>
SymmetricMatrix<TC, VD>::Inverse() const {
	if constexpr(Dimension<=3){ return Comatrix()/Determinant();}
	else { return FromUpperTriangle(this->operator MatrixType().Inverse());}
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
CGSolve(VectorType r) const -> VectorType {
    VectorType x= VectorType::Constant(0);
    const SymmetricMatrix & m = *this;
    VectorType p = r, r0, q;
    ComponentType r2=r.SquaredNorm(), r02=0; //Dummy for r02.
    for(int k=0; k<Dimension && r2>0; ++k){
        if(k>0){
            const ComponentType beta = r2/r02;
            p=r+beta*p;
        }
        q = m*p;
        const ComponentType pq = p.ScalarProduct(q);
        if(pq == 0) return x;
        const ComponentType alpha = r.SquaredNorm()/pq;
        x+=alpha*p;
        r0=r;
        r02=r2;
        r-=alpha*q;
        r2=r.SquaredNorm();
    }
    return x;
}

// ---------------

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Cholevski() const -> CholevskiType {
	std::array<VectorType, Dimension> a;
	for(int i=0; i<Dimension; ++i){
		for(int j=0; j<Dimension; ++j){ // Copy row
			a[i][j] = coef(i,j);}
		for(int j=0; j<i; ++j){
			a[i] -= a[j] * a[i][j]/a[j][j];}
		if(a[i][i]<=0) return {false,MatrixType::Zero()};
		a[i]/=sqrt(a[i][i]);
		}
	return {true,MatrixType::FromRows(a)};
}

// --------------- Norms ------------

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
SquaredFrobeniusNorm() const -> ComponentType {
	ComponentType res=0;
	for(int i=0; i<Dimension; ++i){
		for(int j=0; j<=i; ++j){
			const ComponentType c = coef(i, j);
			if (i==j) {res+=c*c;}
			else {res+=2*c*c;}
		}
	}
	return res;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
FrobeniusNorm() const -> ComponentType {
	using std::sqrt;
	return sqrt(SquaredFrobeniusNorm());
};

/*
template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
SpectralNorm() const -> ComponentType {
	static_assert(Dimension==2);

}*/

#endif /* SymmetricMatrixType_hxx */
