//
//  MatrixType.hxx
//  
//
//  Created by Jean-Marie Mirebeau on 09/02/2018.
//

#ifndef MatrixType_hxx
#define MatrixType_hxx

template<typename TC, size_t VR, size_t VC> template<typename T,typename S>
auto Matrix<TC,VR,VC>::
operator * (const Vector<T, Columns> & u) const -> Vector<S, Rows>
{
	Vector<S, Rows> v;
    v.fill(S(0.));
    for(size_t i=0; i<Rows; ++i){
        for(size_t j=0; j<Columns; ++j)
            v[i]+=this->operator()(i,j)*u[j];
            }
    return v;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Column(int i) const -> OutputVectorType {
	assert(0<=i && i<Columns);
	OutputVectorType result;
	for(int j=0; j<Rows; ++j){result[j] = this->operator()(i,j);}
	return result;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Row(int j) const -> InputVectorType {
	assert(0<=j && j<Rows);
	InputVectorType result;
	for(int i=0; i<Columns; ++i){result[i] = this->operator()(i,j);}
	return result;
}

template<typename TC, size_t VR, size_t VC> template<size_t Columns2> auto Matrix<TC,VR,VC>::
operator * (const Matrix<ComponentType,Columns,Columns2> & m) const
-> Matrix<ComponentType, Rows, Columns2>
{
    Matrix<ComponentType, Rows, Columns2> p;
    p.fill(ComponentType(0));
    for(size_t i=0; i<Rows; ++i)
        for(size_t k=0; k<Columns2; ++k)
            for(size_t j=0; j<Columns; ++j)
                p(i,k) += (this->operator()(i,j))*m(j,k);
    return p;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Transpose() const -> Matrix<ComponentType,Columns,Rows>
{
    Matrix<ComponentType,Columns,Rows> m;
    for(size_t i=0; i<Rows; ++i)
        for(size_t j=0; j<Columns; ++j)
            m(j,i) = this->operator()(i,j);
            return m;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Trace() const -> ComponentType
{
    static_assert(Rows==Columns,"Matrix must be square");
    ComponentType s(0);
    for(size_t i=0; i<Rows; ++i) s+=this->operator()(i,i);
        return s;
}

// ------ Constructors -------
template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Identity() -> Matrix
{
    static_assert(Rows==Columns,"Matrix must be square");
    Matrix m;
    m.fill(ComponentType(0));
    for(size_t i=0; i<Rows; ++i) m(i,i)=1;
    return m;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
FromRows(const std::array<InputVectorType,Rows> & t) -> Matrix
{
    Matrix m;
    for(int i=0; i<Rows; ++i)
        for(int j=0; j<Columns; ++j)
            m(i,j) = t[i][j];
    return m;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
FromColumns(const std::array<OutputVectorType,Columns> & t) -> Matrix
{
    Matrix m;
    for(int i=0; i<Rows; ++i)
        for(int j=0; j<Columns; ++j)
            m(i,j) = t[j][i];
    return m;
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Rotation(ComponentType theta) -> Matrix
{
    static_assert(Rows==2 && Columns==2,"Rotation matrices are two dimensional");
    static_assert( ! std::is_integral<ComponentType>::value,"Rotation matrices have real coefficients");
    const double c=cos(theta), s=sin(theta);
    return FromRows({ InputVectorType{c,-s} , InputVectorType{s,c} });
}

template<typename TC, size_t VR, size_t VC> auto Matrix<TC,VR,VC>::
Random() -> Matrix
{
    Matrix m;
    for(int i=0; i<VR; ++i){
        for(int j=0; j<VC; ++j){
            m(i,j) = -1+2*ComponentType(std::rand())/RAND_MAX;
        }
    }
    return m;
}

template<typename TC, size_t VR, size_t VC> template<typename T> auto Matrix<TC,VR,VC>::
CastCoordinates(const Matrix<T,Rows,Columns> & a0) -> Matrix {
    Matrix a;
    a.data = Vector<ComponentType,Rows*Columns>::CastCoordinates(a0.data);
    return a;
}

template<typename TC, size_t VR, size_t VC>
typename Matrix<TC,VR,VC>::ComponentType
Matrix<TC,VR,VC>::Determinant() const
{
    static_assert(VR==VC,"Matrix must be square");
	const Matrix & m = *this;
	if constexpr(VR==0) {return ComponentType(1.);}
	else if constexpr(VR==1) {return m(0, 0);}
	else if constexpr(VR==2) {return m(0, 0)*m(1, 1) - m(1, 0)*m(0, 1);}
	else if constexpr(VR==3) {
		ComponentType det(0.);
		for (int i = 0; i < 3; ++i) det += m(i, 0)*m((i + 1) % 3, 1)*m((i + 2) % 3, 2) - m(i, 2)*m((i + 1) % 3, 1)*m((i + 2) % 3, 0);
		return det;
	} else {
		// Get largest coefficient, in absolute value
		assert(!std::numeric_limits<TC>::is_integer);
		const size_t d = VR;
		int iMax = 0, jMax = 0;// Dummy initialization
		ComponentType vMax = ComponentType(0);
		for (int i = 0; i < d; ++i)
			for (int j = 0; j < d; ++j) {
				const ComponentType v = std::abs(m(i, j));
				if (v > vMax) {
					iMax = i; jMax = j; vMax = v;
				}
			}
		if (vMax == ComponentType(0)) return ComponentType(0);
		vMax = m(iMax, jMax);
		
		// Pivot
		const size_t dd = (d > 0 ? d - 1 : 0); // Setting dd=d-1 yields recursive template
		Matrix<TC, dd, dd> a; a.fill(TC(0));
		for (int i = 0, ia = 0; i < d; ++i) {
			if (i == iMax) continue;
			const ComponentType delta = m(i, jMax) / vMax;
			for (int j = 0, ja = 0; j < d; ++j) {
				if (j == jMax) continue;
				a(ia, ja) = m(i, j) - delta*m(iMax, j);
				++ja;
			}
			++ia;
		}
		return a.Determinant() *vMax *((iMax + jMax) % 2 == 0 ? 1 : -1);
		
	}
	
}

/*
// -------- Determinant --------

template<typename TC, size_t VR, size_t VC> template<size_t n, typename Dummy>
struct Matrix<TC,VR,VC>::Det {
typedef Matrix<TC,VR,VC> M;
typename M::ComponentType operator()(const M & m){
// Get largest coefficient, in absolute value
assert(!std::numeric_limits<TC>::is_integer);
const size_t d=VR;
int iMax=0,jMax=0;// Dummy initialization
M::ComponentType vMax = TC(0);
for(int i=0; i<d; ++i)
for(int j=0; j<d; ++j){
const ComponentType v = std::abs(m(i,j));
if(v>vMax){
iMax=i; jMax=j; vMax=v;
}
}
if(vMax==TC(0)) return TC(0);
vMax=m(iMax,jMax);

// Pivot
Matrix<TC, d-1, d-1> a; a.fill(TC(0));
for(int i=0, ia=0; i<d; ++i){
if(i==iMax) continue;
const ComponentType delta = m(i,jMax)/vMax;
for(int j=0, ja=0; j<d; ++j){
if(j==jMax) continue;
a(ia,ja)=m(i,j)-delta*m(iMax,j);
++ja;
}
++ia;
}
return a.Determinant() *vMax *((iMax+jMax)%2==0 ? 1 : -1);
}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<0, Dummy> {
typedef Matrix<TC,VR,VC> M;
typename M::ComponentType operator()(const M & m){return TC(1);}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<1, Dummy> {
typedef Matrix<TC,VR,VC> M;
typename M::ComponentType operator()(const M & m){return m(0,0);}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<2, Dummy> {
typedef Matrix<TC,VR,VC> M;
typename M::ComponentType operator()(const M & m){return m(0,0)*m(1,1)-m(1,0)*m(0,1);}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<3, Dummy> {
typedef Matrix<TC,VR,VC> M;
typename M::ComponentType operator()(const M & m){
M::ComponentType det=0;
for(int i=0; i<3; ++i) det+=m(i,0)*m((i+1)%3,1)*m((i+2)%3,2) - m(i,2)*m((i+1)%3,1)*m((i+2)%3,0);
return det;
}
};
*/


// --------- Solve and inverse ---------
template<typename TC, size_t VR, size_t VC>
Matrix<TC,VR,VC>
Matrix<TC,VR,VC>::Inverse() const
{
    static_assert(VR==VC,"Matrix must be square");
    const Matrix & m = *this;
    Matrix b;
	if constexpr(VR==1){
		b(0,0)=ComponentType(1.);
	} else if constexpr(VR==2) {
		b(0,0)= m(1,1);
		b(1,0)=-m(1,0);
		b(0,1)=-m(0,1);
		b(1,1)= m(0,0);
	} else if constexpr(VR==3){
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j)
				b(j,i)=
				m((i+1)%3,(j+1)%3)*m((i+2)%3,(j+2)%3)-
				m((i+1)%3,(j+2)%3)*m((i+2)%3,(j+1)%3);
	} else {
		// Gauss pivot, pasted from solve fct below
            const size_t n=VR;
            Matrix m=*this;
            b=Identity();
            std::array<int,n> i2j, j2i; i2j.fill(-1); j2i.fill(-1);
            for(int j=0; j<n; ++j){
                // Get largest coefficient in column
				ComponentType cMax(0.);
                int iMax=0;
                for(int i=0; i<n; ++i){
                    if(i2j[i]>=0) continue;
                    const ComponentType c = m(i,j);
					using std::abs;
                    if(abs(c)>abs(cMax)){
                        cMax=c; iMax=i;}
                }
                i2j[iMax]=j;
                j2i[j]=iMax;
                assert(cMax!=0); // Otherwise, matrix is not invertible
                
                // Remove line from other lines, while performing likewise on b
                for(int i=0; i<n; ++i){
                    if(i2j[i]>=0) continue;
                    const ComponentType r = m(i,j)/cMax;
                    for(int k=j+1; k<n; ++k){
                        m(i,k)-=m(iMax,k)*r;}
                    for(int l=0; l<VR; ++l){
                        b(i,l)-=b(iMax,l)*r;}
                }
            }
            // Solve remaining triangular system
            Matrix a;
            for(int j=n-1; j>=0; --j){
                const int i=j2i[j];
                //ComponentType & r = a[j];
                for(int l=0; l<VR; ++l){
                a(j,l)=b(i,l);
                for(int k=j+1; k<n; ++k){
                    a(j,l)-=a(k,l)*m(i,k);}
                a(j,l)/=m(i,j);
                }
            }
            return a;
	}  // if dimension
	
    const ComponentType d=Determinant();
	if constexpr(std::is_same_v<decltype(d!=ComponentType(0.)),bool>) assert(d!=0); // exception ?
    return b/d;
}

template<typename TC, size_t VR, size_t VC> auto
Matrix<TC,VR,VC>::Solve(OutputVectorType b) const -> InputVectorType {
    // A basic Gauss pivot
    static_assert(VR==VC,"Matrix must be square");
    const size_t n=VR;
    std::array<int,n> i2j, j2i; i2j.fill(-1); j2i.fill(-1);
    Matrix m = this;
    for(int j=0; j<n; ++j){
        // Get largest coefficient in column
        ComponentType cMax=0;
        int iMax=0;
        for(int i=0; i<n; ++i){
            if(i2j[i]>=0) continue;
            const ComponentType c = m(i,j);
            if(std::abs(c)>std::abs(cMax)){
                cMax=c; iMax=i;}
        }
        i2j[iMax]=j;
        j2i[j]=iMax;
        assert(cMax!=0); // Matrix is not invertible
        
        // Remove line from other lines, while performing likewise on b
        for(int i=0; i<n; ++i){
            if(i2j[i]>=0) continue;
            const ComponentType r = m(i,j)/cMax;
            for(int k=j+1; k<n; ++k){
                m(i,k)-=m(iMax,k)*r;}
            b[i]-=b[iMax]*r;
        }
    }
    // Solve remaining triangular system
    InputVectorType a;
    for(int j=n-1; j>=0; --j){
        const int i=j2i[j];
        ComponentType & r = a[j];
        r=b[i];
        for(int k=j+1; k<n; ++k){
            r-=a[k]*m(i,k);}
        r/=m(i,j);
    }
    return a;
}

#endif /* MatrixType_hxx */
