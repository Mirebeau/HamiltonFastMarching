//
//  ScalarTraits.h
//  SimdTest
//
//  Created by Jean-Marie Mirebeau on 12/04/2019.
//  Copyright © 2019 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef ScalarTraits_h
#define ScalarTraits_h

/*
 This header implements basic features for working with compound types in place of scalars.
 It was used with "xsimd.h" library.
 */

#include <type_traits>

template<typename T> constexpr bool is_simple_v
= std::is_same_v<decltype(std::declval<T>()==std::declval<T>()),bool>;

/*
template<typename T, bool b = is_simple_v<T> > struct get_type;
template<typename T> struct get_type <T,true> {using type=T;};
template<typename T> struct get_type <T,false> {using type= typename T::type;};
template<typename T> using get_type_t = typename get_type<T>::type;
 */

#endif /* ScalarTraits_h */
