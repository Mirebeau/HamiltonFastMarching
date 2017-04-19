// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef PeriodicGrid_hxx
#define PeriodicGrid_hxx

// Basic conversions
template<typename T> auto PeriodicGrid<T>::
LinearFromIndex(IndexCRef index) const -> DiscreteType {
    return arr.Convert(index);}

template<typename T> auto PeriodicGrid<T>::
IndexFromLinear(DiscreteType linearIndex) const -> IndexType {
    return arr.Convert(linearIndex);}


template<typename T> auto PeriodicGrid<T>::
IndexFromPoint(const PointType & p) const -> IndexType {
    IndexType result;
    for(int i=0; i<Dimension; ++i){
        result[i] = std::floor(p[i]);}
    return result;
}

template<typename T> auto PeriodicGrid<T>::
PointFromIndex(const IndexType & p) const -> PointType {
    PointType result;
    for(int i=0; i<Dimension; ++i){
        result[i]=p[i]+0.5;}
    return result;
}

// Boundary conditions
template<typename T> constexpr bool PeriodicGrid<T>::
MayReverse(DiscreteType i) {
//    assert(0<=i && i<Dimension); // Some compilers complain about assert in constexpr
    return
    Traits::boundaryConditions[i]==Boundary::Sphere2_0 ||
    Traits::boundaryConditions[i]==Boundary::Sphere3_0;
}

template<typename T> PeriodicGrid<T>::
PeriodicGrid(IndexCRef dims){
    arr.dims=dims;
    const auto & bc=Traits::boundaryConditions;
    for(int i=0; i<Dimension; ++i){
        const DiscreteType n=dims[i];
        if(n<=0)
            ExceptionMacro("PeriodicGrid error : non-positive dimension");
        
        switch(bc[i]){
            case Boundary::Sphere2_0:{
                if(i>Dimension-2) ExceptionMacro("PeriodicGrid error : Sphere2_0 must be followed by Sphere2_1");
                if(bc[i+1]==Boundary::Sphere2_1){
                    if(dims[i+1]!=2*n && dims[i+1]!=4*n) ExceptionMacro("PeriodicGrid error : inconsistent Sphere2 dimensions");
                } else if(bc[i+1]==Boundary::Sphere2_0) {
                    if(dims[i+1]!=n && dims[i+1]!=2*n) ExceptionMacro("PeriodicGrid error : inconsistent Sphere2 dimensions");
                } else {ExceptionMacro("PeriodicGrid error : Sphere2_0 must be followed by Sphere2_1");}
                break;}
                
            case Boundary::Sphere3_0:{
                if(i>Dimension-3) ExceptionMacro("PeriodicGridError: Sphere3_0 must be followed by Sphere3_1|Proj3_1, Sphere3_2");
                if(bc[i+1]==Boundary::Sphere3_1 && bc[i+2]==Boundary::Sphere3_2){
                    if(dims[i+1]!=4*n || dims[i+2]!=4*n) ExceptionMacro("PeriodicGrid error : inconsistent Sphere3 dimensions");
                } else if (bc[i+1]==Boundary::Proj3_1 && bc[i+2]==Boundary::Sphere3_2){
                    if(dims[i+1]!=2*n || dims[i+2]!=4*n) ExceptionMacro("PeriodicGrid error : inconsistent Sphere3 dimensions");
                } else {ExceptionMacro("PeriodicGridError: Sphere3_0 must be followed by Sphere3_1|Proj3_1, Sphere3_2");}
                break;}
                
            case Boundary::Sphere2_1:{
                if(i<=0 || bc[i-1]!=Boundary::Sphere2_0)
                    ExceptionMacro("PeriodicGrid error : Sphere2_1 must be preceded by Sphere2_0");
                break;}
            case Boundary::Sphere3_1:{
                if(i<=0 || bc[i-1]!=Boundary::Sphere3_0)
                    ExceptionMacro("PeriodicGrid error : Sphere3_1 must be preceded by Sphere3_0");
                break;}
            case Boundary::Proj3_1:{
                if(i<=0 || bc[i-1]!=Boundary::Sphere3_0)
                    ExceptionMacro("PeriodicGrid error : Proj3_1 must be preceded by Sphere3_0");
                break;}
            case Boundary::Sphere3_2:{
                if(i<=1 || bc[i-2]!=Boundary::Sphere3_0)
                    ExceptionMacro("PeriodicGrid error : Sphere3_2 must be preceded by Sphere3_0");
                break;}
                
            case Boundary::Closed:
            case Boundary::Periodic:
                break;
        }
    }
}


// Periodization
template<typename T> auto PeriodicGrid<T>::
Periodize(IndexType & point) const -> ReverseFlag {
    ReverseFlag flipped;
    const auto & dims = arr.dims;
    for(int i=0; i<Dimension; ++i){
        switch(Traits::boundaryConditions[i]) {
                
            case Boundary::Closed:
                if(point[i]<0 || point[i]>=dims[i]) {
                    flipped[Dimension]=1;
                    return flipped;
                }
                break;
                
            case Boundary::Periodic:
            case Boundary::Sphere2_1:
            case Boundary::Sphere3_1:
            case Boundary::Sphere3_2:
                point[i] = PosMod(point[i], dims[i]);
                break;
                
            case Boundary::Sphere2_0:
                point[i]=PosMod(point[i], 2*dims[i]);
                if(point[i]>=dims[i]){
                    point[i] = 2*dims[i]-point[i]-1;
                    flipped[i]=1;
                    if(T::boundaryConditions[i+1]==Boundary::Sphere2_0) {
                        // Hack to handle arbitrary dimensional spheres,
                        // but this parametrization is ill conditioned in dim n>2.
                        point[i+1]=point[i+1]+dims[i+1];
                    } else if(T::boundaryConditions[i+1]==Boundary::Sphere2_1) {
                        point[i+1]=point[i+1]+dims[i+1]/2;
                    } else {assert(false);}
                }
                assert(0<=point[i] && point[i]<dims[i]);
                break;
                
            case Boundary::Sphere3_0:{
                const DiscreteType n=dims[i]; //equiv pi/2
                point[i]=PosMod(point[i], 4*n);
                if(point[i]>=2*n){
                    point[i]-=2*n;
                    point[i+1]+=2*n;
                    point[i+2]+=2*n;
                }
                if(point[i]>=n){
                    point[i]=2*n-point[i]-1;
                    point[i+1]+=2*n;
                }
            }
                
            case Boundary::Proj3_1:{
                const DiscreteType n=dims[i]; //equiv pi
                point[i]=PosMod(point[i], 2*n);
                if(point[i]>=n){
                    point[i]-=n;
                    point[i+1]+=n;
                }
            }
        }
    }
    return flipped;
}

template<typename T> auto PeriodicGrid<T>::
Periodize(PointType & point) const -> ReverseFlag {
    ReverseFlag flipped;
    const auto & dims = arr.dims;
    for(int i=0; i<Dimension; ++i){
        switch (T::boundaryConditions[i]) {
            case Boundary::Closed:
                if(point[i]<0 || point[i]>=dims[i]) {
                    flipped[Dimension]=1;
                    return flipped;
                }
                break;
                
            case Boundary::Periodic:
            case Boundary::Sphere2_1:
            case Boundary::Sphere3_1:
            case Boundary::Sphere3_2:
                point[i] = fPosMod(point[i],(ScalarType)dims[i]);
                break;
                
            case Boundary::Sphere2_0:
                assert(i<Dimension-1);
                point[i] = fPosMod(point[i],(ScalarType)2*dims[i]);
                if(point[i]>=dims[i]){
                    point[i] = 2*dims[i]-point[i]; // !! Note the missing -1 w.r.t discrete
                    flipped[i]=1;
                    if(T::boundaryConditions[i+1]==Boundary::Sphere2_0) {
                        point[i+1]=point[i+1]+dims[i+1];
                    } else if(T::boundaryConditions[i+1]==Boundary::Sphere2_1) {
                        assert(dims[i+1]%2==0);
                        point[i+1]=point[i+1]+dims[i+1]/2;
                    } else {assert(false);}
                }
                assert(0<=point[i] && point[i]<dims[i]);
                break;
                
            case Boundary::Sphere3_0:{
                const DiscreteType n=dims[i]; //equiv pi/2
                point[i]=fPosMod(point[i], (ScalarType)4*n);
                if(point[i]>=2*n){
                    point[i]-=2*n;
                    point[i+1]+=2*n;
                    point[i+2]+=2*n;
                }
                if(point[i]>=n){
                    point[i]=2*n-point[i];
                    point[i+1]+=2*n;
                }
            }
                
            case Boundary::Proj3_1:{
                const DiscreteType n=dims[i]; //equiv pi
                point[i]=fPosMod(point[i], (ScalarType)2*n);
                if(point[i]>=n){
                    point[i]-=n;
                    point[i+1]+=n;
                }
            }
                
        }
    }
    return flipped;
}


#endif /* PeriodicGrid_hxx */
