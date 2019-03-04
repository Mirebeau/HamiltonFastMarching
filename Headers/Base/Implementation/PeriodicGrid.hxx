// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef PeriodicGrid_hxx
#define PeriodicGrid_hxx

// ------- Methods of child class transform -------

// Boundary conditions
template<typename T> constexpr bool PeriodicGrid<T>::Transform::
MayReverse(DiscreteType i) {
//    assert(0<=i && i<Dimension); // Some compilers complain about assert in constexpr
    return
    Traits::boundaryConditions[i]==Boundary::Sphere2_0 ||
    Traits::boundaryConditions[i]==Boundary::Sphere3_0;
}

template<typename T> template<typename TVec> void PeriodicGrid<T>::
Transform::PullVector(TVec & v) const {
    for(DiscreteType i=0; i<Dimension; ++i){
        if(MayReverse(i) && reverseFlag[i])
            v[i]*=-1;
    }
}

// ----- Neighbors on the grid -----
template<typename T> auto PeriodicGrid<T>::
Neighbors(const PointType & p, bool requireTrivial) const -> NeighborsType {
	auto result = Superclass::Neighbors(p);
	ScalarType wSum = 0.;
	for(auto & [q,w] : result){
		const Transform transform = PeriodizeNoBase(q);
		if(transform.IsTrivial() || (transform.IsTrivial() && !requireTrivial) ) {wSum+=w;}
		else {w=0.;}
	}
	if(wSum==0)
		ExceptionMacro("Error : point " << p << " has no valid neighbors in the grid domain");
	return result;
}


/*
 //	template<typename T, typename F> T Interpolate(const F &, bool trivialOnly) const;

 template<typename TTraits> template<typename T, typename F> auto PeriodicGrid<TTraits>::
Interpolate(const F &, bool TrivialOnly) const {
	const auto

}*/

// --------- Main methods -------


template<typename T> PeriodicGrid<T>::
PeriodicGrid(IndexCRef dims): Superclass(dims){
    const std::string err = "PeriodicGrid error : ";
    const auto & bc=Traits::boundaryConditions;
    for(int i=0; i<Dimension; ++i){
        const DiscreteType n=dims[i];
        
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
                if(dims[i]%2 != 0) ExceptionMacro(err << "Hi there");
                break;}
            case Boundary::Sphere2_Hopf:{
                if(i<=1 || bc[i-1]!=Boundary::Sphere2_1)
                    ExceptionMacro("PeriodicGrid error : Sphere2_Hopf must be preceded by Sphere2_1");
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
Periodize(IndexType & point, IndexCRef) const -> Transform {
    return PeriodizeNoBase(point);}

template<typename T> auto PeriodicGrid<T>::
Periodize(PointType & point, const PointType &) const -> Transform {
    return PeriodizeNoBase(point);}

template<typename T> auto PeriodicGrid<T>::
PeriodizeNoBase(IndexType & point) const -> Transform {
    Transform flipped;
    const auto & dims = this->arr.dims;
    for(int i=0; i<Dimension; ++i){
        switch(Traits::boundaryConditions[i]) {
                
            case Boundary::Closed:
                if(point[i]<0 || point[i]>=dims[i]){
                    flipped.Invalidate();
                    return flipped;
                }
                break;
                
            case Boundary::Periodic:
            case Boundary::Sphere2_1:
            case Boundary::Sphere3_1:
            case Boundary::Sphere3_2:
            case Boundary::Sphere2_Hopf:
                point[i] = PosMod(point[i], dims[i]);
                break;
                
            case Boundary::Sphere2_0:
                point[i]=PosMod(point[i], 2*dims[i]);
                if(point[i]>=dims[i]){
                    point[i] = 2*dims[i]-point[i]-1;
                    flipped.reverseFlag[i]=1;
                    if(T::boundaryConditions[i+1]==Boundary::Sphere2_0) {
                        // Hack to handle arbitrary dimensional spheres,
                        // but this parametrization is ill conditioned in dim n>2.
                        point[i+1]=point[i+1]+dims[i+1];
                    } else if(T::boundaryConditions[i+1]==Boundary::Sphere2_1) {
                        point[i+1]=point[i+1]+dims[i+1]/2;
                        if(i<=Dimension-3 && T::boundaryConditions[i+2]==Boundary::Sphere2_Hopf){
                            point[i+2]=point[i+2]+dims[i+2]/2;
                        }
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
                    flipped.reverseFlag[i]=1;
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
PeriodizeNoBase(PointType & point) const -> Transform {
    Transform flipped;
    const auto & dims = this->arr.dims;
    for(int i=0; i<Dimension; ++i){
        switch (T::boundaryConditions[i]) {
            case Boundary::Closed:
                if(point[i]<0 || point[i]>=dims[i]) {
                    flipped.Invalidate();
                    return flipped;
                }
                break;
                
            case Boundary::Periodic:
            case Boundary::Sphere2_1:
            case Boundary::Sphere3_1:
            case Boundary::Sphere3_2:
            case Boundary::Sphere2_Hopf:
                point[i] = fPosMod(point[i],(ScalarType)dims[i]);
                break;
                
            case Boundary::Sphere2_0:
                assert(i<Dimension-1);
                point[i] = fPosMod(point[i],(ScalarType)2*dims[i]);
                if(point[i]>=dims[i]){
                    point[i] = 2*dims[i]-point[i]; // !! Note the missing -1 w.r.t discrete
                    flipped.reverseFlag[i]=1;
                    if(T::boundaryConditions[i+1]==Boundary::Sphere2_0) {
                        point[i+1]=point[i+1]+dims[i+1];
                    } else if(T::boundaryConditions[i+1]==Boundary::Sphere2_1) {
                        assert(dims[i+1]%2==0);
                        point[i+1]=point[i+1]+dims[i+1]/2;
                        if(i<=Dimension-3 && T::boundaryConditions[i+2]==Boundary::Sphere2_Hopf){
                            point[i+2]=point[i+2]+dims[i+2]/2;
                        }
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
                    flipped.reverseFlag[i]=1;
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
