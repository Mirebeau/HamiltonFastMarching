//
//  AlignedBillard.hxx
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 15/08/2018.
//

#ifndef AlignedBillard_hxx
#define AlignedBillard_hxx




// ----  Periodize ------
template<typename T> auto
AlignedBillardGrid<T>::Periodize(IndexType & target, IndexCRef base)
const -> Transform {
    if(!inside.empty() && inside.InRange(target) && inside(target)) {
        return Transform();}
    PointType targetP=this->PointFromIndex(target), baseP = this->PointFromIndex(base);
    const Transform & result = Periodize(targetP,baseP);
    target = this->IndexFromPoint(targetP);
    assert(this->arr.InRange(target));
    return result;
}

template<typename T> auto
AlignedBillardGrid<T>::Periodize(PointType & target, PointType base)
const -> Transform {
    assert(target!=base);
    assert(std::fabs(LoopIndex(base)-1)<0.01);
    
    const VectorType v = target-base;
    ScalarType s=0;
    while(s<1){
        s=Cut(base,v,s);}
    
    target = base+v;
    assert(std::fabs(LoopIndex(target)-1.)<0.01);
    
    return Transform();
}


template<typename T> auto
AlignedBillardGrid<T>::PeriodizeNoBase(PointType & target)
const -> Transform {
    if(std::fabs(LoopIndex(target)-1.)<0.01){return Transform();}
    Transform result; result.Invalidate(); return result;
}


template<typename T> auto
AlignedBillardGrid<T>::PeriodizeNoBase(IndexType & target)
const -> Transform {
    assert(!inside.empty());
    if(inside.InRange(target) && inside(target)) {return Transform();}
    Transform result; result.Invalidate(); return result;
}
    
// -------- Geometry -------

template<typename T> auto
AlignedBillardGrid<T>::Cut(PointType & p, const VectorType & v, ScalarType s0)
const -> ScalarType {
    VectorType trans = VectorType::Constant(0.);
    ScalarType s = 1;
    for(const EdgeType & e : edges){
        const ScalarType sTest = Cut(p,v,e);
        if(sTest < s0 || sTest >= s) continue;
        s=sTest;
        trans = e.trans;
    }
    p+=trans;
    return s;
}

template<typename T> auto
AlignedBillardGrid<T>::Cut(const PointType & p, const VectorType & v, const EdgeType & e)
const -> ScalarType {
    const ScalarType delta = LinearAlgebra::Determinant(v, e.v);
    if(delta<=0) return Traits::Infinity();
    const VectorType w = p-e.p;
    const ScalarType t = LinearAlgebra::Determinant(v,w)/delta;
    if(!(0<=t && t<=1)) return Traits::Infinity();
    return LinearAlgebra::Determinant(e.v,w)/delta;
}

// ------ Points in domain ------

template<typename T> auto
AlignedBillardGrid<T>::LoopIndex(const PointType & p) const -> ScalarType {
    ScalarType theta = 0;
    for(const EdgeType & e : edges){
        const VectorType v0 = e.p - p;
        const VectorType v1 = v0 + e.v;
        theta += atan2(LinearAlgebra::Determinant(v0, v1),
                       v0.ScalarProduct(v1));
    }
    return theta/(2.*Traits::mathPi);
}


template<typename T> void
AlignedBillardGrid<T>::SetPeriodizeLazy() const {
    inside.dims = this->arr.dims;
    inside.resize(inside.dims.Product());
    for(DiscreteType i=0; i<inside.size(); ++i){
        const PointType p =this->PointFromIndex(this->IndexFromLinear(i));
        const ScalarType loop = LoopIndex(p);
        if(std::fabs(loop)<0.01){inside[i]=false;}
        else if(std::fabs(1-loop)<0.01){inside[i]=true;}
        else {ExceptionMacro("AlignedBillardGrid error: index of point "
                             << p << " is not in {0,1}.");}
    }
}


#endif /* AlignedBillard_hxx */
