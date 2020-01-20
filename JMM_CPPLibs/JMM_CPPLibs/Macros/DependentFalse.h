//
//  DependentFalse.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 12/12/2019.
//

#ifndef DependentFalse_h
#define DependentFalse_h

template<typename T> struct dependent_false : std::false_type {};
template<typename T> constexpr bool dependent_false_v = dependent_false<T>::value;

#endif /* DependentFalse_h */
