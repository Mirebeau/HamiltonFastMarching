//
//  VectorType.hxx
//  
//
//  Created by Jean-Marie Mirebeau on 12/02/2018.
//

#ifndef VectorType_hxx
#define VectorType_hxx

template<typename TC, size_t VD> auto Vector<TC,VD>::
Normalized() const -> Vector {
    ComponentType norm = Norm();
    return norm==0 ? *this : (*this)/norm;
}

template<typename TC, size_t VD> auto Vector<TC,VD>::
RandomUnit() -> Vector {
    assert( !std::is_integral<ComponentType>::value );
    Vector v;
    do {
        // Draw a random unit vector uniformly in the cube.
        for(int i=0; i<Dimension; ++i){
            v[i] = -1+2*double(rand())/RAND_MAX;}
        // Reject if not in the ball
    } while(v.SquaredNorm()>1);
    return v.Normalized();
}

#endif /* VectorType_hxx */
