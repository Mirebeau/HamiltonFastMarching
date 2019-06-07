//
//  MapWeightedSum.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/02/2019.
//

#ifndef MapWeightedSum_h
#define MapWeightedSum_h

/*
 - f is a functor
 - data is an array of pairs of arguments to f and weights
 
 returns sum_{ [x,w] : data} w*f(x)
 
 * Omits evaluation of f(x) when w=0. (Allows possibly invalid arguments with zero weight.)
 * Assumes first element of data is valid argument.
 
 */

// -> decltype(data.begin()->second * f(data.begin()->first) )
template<typename T, typename F, typename ArgWeightData>
T MapWeightedSum(const F & f, const ArgWeightData & data)
{
	auto it = data.begin(), end = data.end();
	assert(it!=end);
	T result = it->second * f(it->first);
	for(++it; it!=end; ++it){
		if(it->second!=0){
			result += it->second * f(it->first);}
	}
	return result;
}

#endif /* MapWeightedSum_h */
