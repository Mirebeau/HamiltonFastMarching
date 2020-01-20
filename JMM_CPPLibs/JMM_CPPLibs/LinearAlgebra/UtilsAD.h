//
//  UtilsAD.h
//  TestTTINorm
//
//  Created by Jean-Marie Mirebeau on 17/01/2020.
//  Copyright © 2020 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef UtilsAD_h
#define UtilsAD_h

namespace LinearAlgebra {

template<typename TScalar> struct UtilsAD {
	using ScalarType = TScalar;
	using ADType = ScalarType;
	static_assert(std::numeric_limits<ScalarType>::is_specialized, "Unsupported type");
	static ScalarType Scalar(const ScalarType & t){return t;}
};

template<typename T> typename UtilsAD<T>::ScalarType
RemoveAD(const T & t){return UtilsAD<T>::Scalar(t);}

} // Namespace LinearAlgebra
#endif /* UtilsAD_h */
