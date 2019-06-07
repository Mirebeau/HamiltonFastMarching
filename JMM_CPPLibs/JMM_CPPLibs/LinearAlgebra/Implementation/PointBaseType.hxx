//
//  PointBaseType.hxx
//  
//
//  Created by Jean-Marie Mirebeau on 09/02/2018.
//

#ifndef PointBaseType_hxx
#define PointBaseType_hxx

/* // Makes things ambiguous
 explicit PointBase(std::initializer_list<ComponentType> data){
 assert(data.size()==Dimension);
 std::copy(data.begin(),data.end(),this->begin());
 };*/

template<typename TC, size_t VD>
bool
PointBase<TC,VD>::IsNonNegative() const {
    for(int i=0; i<Dimension; ++i)
        if(this->operator[](i)<0) return false;
    return true;
}

template<typename TC, size_t VD>
bool
PointBase<TC,VD>::IsPositive() const {
    for(int i=0; i<Dimension; ++i)
        if(this->operator[](i)<=0) return false;
    return true;
}

template<typename TC, size_t VD>
bool
PointBase<TC,VD>::IsFinite() const {
    for(int i=0; i<Dimension; ++i)
        if(!(std::fabs(this->operator[](i))<std::numeric_limits<ComponentType>::infinity()))
            return false;
    return true;
}

template<typename TC, size_t VD>
typename PointBase<TC,VD>::ComponentType
PointBase<TC,VD>::Sum() const {
    ComponentType sum(0);
    for(int i=0; i<Dimension; ++i)
        sum += this->operator[](i);
    return sum;
}

template<typename TC, size_t VD>
typename PointBase<TC,VD>::ComponentType
PointBase<TC,VD>::Product() const {
    ComponentType prod(1);
    for(int i=0; i<Dimension; ++i)
        prod *= this->operator[](i);
    return prod;
}

template<typename TC, size_t VD>
struct PointBase<TC, VD>::LexicographicCompare {
    bool operator()(const PointBase & p, const PointBase & q) const {
        for(auto pi=p.begin(), qi=q.begin(); pi!=p.end(); ++pi, ++qi){
            if(*pi < *qi) return true;
            if(*pi > *qi) return false;
        }
        return false; // strict ordering
    }
};

#endif /* PointBaseType_hxx */
