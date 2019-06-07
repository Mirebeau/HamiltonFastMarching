//
//  GetComponent.h
//  TestFromComponentsIterator
//
//  Created by Jean-Marie Mirebeau on 22/09/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef GetComponent_h
#define GetComponent_h

// Access to the components of an std::array, scalar (seen as array of length 1), or std::pair of these
template<typename T, typename C>
struct GetComponent {
	static constexpr size_t size() {return T().size();}
	static const C & Get(const T & t, size_t i) {return t[i];}
	static C & Get(T & t, size_t i) {return t[i];}
};

template<typename TA,typename TB, typename C>
struct GetComponent<std::pair<TA,TB>, C> {
	typedef std::pair<TA,TB> T;
	static constexpr size_t size() {return sizea()+GetComponent<TB,C>::size();}
	static const C & Get(const T & t, size_t i) {
		return i<sizea() ? GetComponent<TA,C>::Get(t.first,i) : GetComponent<TB,C>::Get(t.second,i-sizea());}
	static C & Get(T & t, size_t i) {
		return i<sizea() ? GetComponent<TA,C>::Get(t.first,i) : GetComponent<TB,C>::Get(t.second,i-sizea());}
protected:
	static constexpr size_t sizea() {return GetComponent<TA,C>::size();}
};

template<typename T>
struct GetComponent<T, T> {
	static constexpr size_t size() {return 1;}
	static const T & Get(const T & t, size_t i) {assert(i==0); return t;}
	static T & Get(T & t, size_t i) {assert(i==0); return t;}
};


#endif /* GetComponent_h */
