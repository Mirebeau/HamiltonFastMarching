//
//  FromComponentsIterator.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 21/09/2018.
//

#ifndef FromComponentsIterator_h
#define FromComponentsIterator_h
#include <limits>
#include "../LinearAlgebra/FriendOperators.h"
#include "../Macros/RedeclareTypes.h"
#include "GetComponent.h"

// Transform an iterator to the components of an object type into an iterator to these objects.
// Acts (mostly) like a pointer if the types are identical
// Usage ex : FromComponentsIterator<std::array<double,3>, double>, or const double for read only
template<typename TT, typename TComp>
struct FromComponentsIterator
: LinearAlgebra::offsettable<FromComponentsIterator<TT,TComp>, long> {
	typedef TT T;
	typedef typename std::remove_const<TComp>::type ComponentType;
	typedef TComp ElementType;
	typedef FromComponentsIterator FCI;
	typedef GetComponent<T,ComponentType> GetComp;
	
	// Constructor, and factory from const pointer.
	FromComponentsIterator(ElementType * it_):it(it_){};
	static constexpr size_t nComp() {return GetComp::size();}

	// Access to pointee
	const T operator * () const {
		T t{};
		for(size_t i=0; i<nComp(); ++i){GetComp::Get(t,i)=it[i];}
		return t;
	}
	void Set(const T & t){
		for(size_t i=0; i<nComp(); ++i){it[i]=GetComp::Get(t,i);}
	}
	
	// Increment, decrement
	FCI & operator ++(){it+=nComp(); return self();}
	FCI operator ++(int){FCI r=self(); ++self(); return r;}
	FCI & operator --(){it-=nComp(); return self();}
	FCI operator --(int){FCI r=self(); --self(); return r;}
	
	// Offset by long
	FCI & operator +=(long i){it+=nComp()*i; return self();}
	FCI & operator -=(long i){return self()+=-i;}
	const T operator[](long i) const {return *(self()+i);}
protected:
	ElementType * it = nullptr;
	FCI & self(){return *this;}
	const FCI & self() const {return *this;}
};

#endif /* FromComponentsIterator_h */
