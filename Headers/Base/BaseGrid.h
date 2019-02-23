//
//  BaseGrid.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 15/08/2018.
//

#ifndef BaseGrid_h
#define BaseGrid_h

#include "JMM_CPPLibs/LinearAlgebra/ArrayType.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorType.h"

template<int VDim, typename TScalar, typename TDiscrete, typename TShort> struct BaseGrid {
    static const int Dimension = VDim;
    typedef TScalar ScalarType;
    typedef TDiscrete DiscreteType;
    typedef TShort ShortType;
    
    typedef LinearAlgebra::Point<ScalarType, Dimension> PointType;
    typedef LinearAlgebra::Vector<ScalarType, Dimension> VectorType;
    
    /// Multi-index of a grid point
    typedef LinearAlgebra::Point<DiscreteType, Dimension> IndexType;
    typedef LinearAlgebra::Vector<DiscreteType, Dimension> IndexDiff;
    /// Compact storage of offsets between neighbor indices
    typedef LinearAlgebra::Vector<ShortType, Dimension> OffsetType;
    
    typedef const IndexType & IndexCRef;
    template<typename T, size_t n> using Array=LinearAlgebra::Array<T,n>;
    
    DiscreteType LinearFromIndex(IndexCRef index) const {return arr.Convert(index);}
    IndexType    IndexFromLinear(DiscreteType linearIndex) const {return arr.Convert(linearIndex);}
    static PointType    PointFromIndex(IndexCRef);
    static IndexType    IndexFromPoint(const PointType &);

    BaseGrid(IndexCRef dims) {arr.dims = dims;
        if(!dims.IsPositive()) ExceptionMacro("BaseGrid error : domain size must be positive");}
	
	using NeighborsType = std::array<std::pair< IndexType,ScalarType>, 1<<Dimension >;
	static NeighborsType Neighbors(const PointType &); // Neighbors on the cartesian grid Z^d
protected:
    Array<ScalarType, Dimension> arr;
};


// A second class is used to reparametrize the domain
template<typename TPoint, typename TVec> struct ParamInterface_ {
    typedef TPoint PointType;
    typedef TVec VectorType;
    virtual PointType   ADim(   const PointType & p)    const=0;
    virtual VectorType  ADim(   const VectorType & v)   const = 0;
    virtual PointType   ReDim(  const PointType & p)    const=0;
    virtual VectorType  ReDim(  const VectorType & v)   const=0;
    virtual ~ParamInterface_(){};
};

// -------------------- Implementation ------------------

// Basic conversions
template<int VD, typename TS, typename TD, typename TSh> auto BaseGrid<VD,TS,TD,TSh>::
IndexFromPoint(const PointType & p) -> IndexType {
    IndexType result;
    for(int i=0; i<Dimension; ++i){
        result[i] = (DiscreteType)std::floor(p[i]);}
    return result;
}

template<int VD, typename TS, typename TD, typename TSh> auto BaseGrid<VD,TS,TD,TSh>::
PointFromIndex(const IndexType & p) -> PointType {
    PointType result;
    for(int i=0; i<Dimension; ++i){
        result[i]=p[i]+0.5;}
    return result;
}

// Neighbors

template<int VD, typename TS, typename TD, typename TSh> auto BaseGrid<VD,TS,TD,TSh>::
Neighbors(const PointType & p) -> NeighborsType {
	const IndexType nearest = IndexFromPoint(p);
	IndexType direction;
	PointType weight;
	for(int i=0; i<Dimension; ++i){
		const ScalarType diff = p[i]-nearest[i]-0.5;
		direction[i] = diff>0 ? 1 : -1;
		weight[i] = std::abs(diff);
	}
	
	NeighborsType result;
	for(int i=0; i< (1<<Dimension); ++i){
		IndexType p = nearest;
		ScalarType w = 1.;
		for(int j=0; j<Dimension; ++j){
			const bool b = (bool)( (i>>j) & 1 );
			if(b){
				w *= weight[j];
				p[j]+=direction[j];
			} else {
				w*= 1.-weight[j];
			}
		}
		result[i] = {p,w};
	}
	return result;
}

#endif /* BaseGrid_h */
