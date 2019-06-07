//
//  MultiVector.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 22/02/2019.
//

#ifndef MultiVector_h
#define MultiVector_h

#include <vector>
#include "RangeAccessor.h"

template<typename TValue, typename TSize=long>
struct MultiVector {
	using ValueType = TValue;
	using SizeType = TSize;
	using ContainerType = std::vector<ValueType>;
	using Iterator = typename ContainerType::iterator;
	using ConstIterator = typename ContainerType::const_iterator;
	using RangeType = RangeAccessor<Iterator>;
	using ConstRangeType = RangeAccessor<ConstIterator>;
	using KeyValuePair = std::pair<SizeType,ValueType>;

	ContainerType values;
	std::vector<SizeType> splits;
	
	RangeType operator[](SizeType);
	ConstRangeType operator[](SizeType) const;
	
	RangeType front(){assert(!empty()); return operator[](0);}
	ConstRangeType front() const {assert(!empty()); return operator[](0);}
	RangeType back() {assert(!empty()); return operator[](size()-1);}
	ConstRangeType back() const {assert(!empty()); return operator[](size()-1);}

	MultiVector(){splits.push_back(0);}

	SizeType size() const;
	void clear();
	bool empty(){return splits.size()==1;}
	void pop_back();
	template<typename It> void push_back(const RangeAccessor<It> &);
	void insert(const std::vector<KeyValuePair> &, SizeType=0); // must be sorted
protected:
	bool IsValid(SizeType) const;
};

template<typename TV,typename TS> bool
MultiVector<TV,TS>::IsValid(SizeType index) const {
	return (index+1<splits.size() &&
	splits[index]<=splits[index+1] &&
	splits[index+1]<=values.size());
}

template<typename TV,typename TS> auto
MultiVector<TV,TS>::operator[](SizeType index)
-> RangeAccessor<Iterator> {
	assert(IsValid(index));
	return {values.begin()+splits[index], values.begin()+splits[index+1]};
}

template<typename TV,typename TS> auto
MultiVector<TV,TS>::operator[](SizeType index) const
-> RangeAccessor<ConstIterator> {
	assert(IsValid(index));
	return {values.begin()+splits[index], values.begin()+splits[index+1]};
}

template<typename TV,typename TS> auto
MultiVector<TV,TS>::size() const -> SizeType {
	assert(!splits.empty());
	return splits.size()-1;
}

template<typename TV,typename TS> void
MultiVector<TV,TS>::pop_back(){
	assert(size()>0);
	splits.pop_back();
	values.resize(splits.back());
}

template<typename TV,typename TS> void
MultiVector<TV,TS>::clear(){
	values.clear();
	splits.resize(1);
	assert(splits[0]==0);
}

template<typename TV,typename TS> template<typename It> void
MultiVector<TV,TS>::push_back(const RangeAccessor<It> & data){
	values.reserve(values.size()+data.size());
	for(const ValueType & val : data){values.push_back(val);}
	splits.push_back((SizeType)values.size());
}

template<typename TV,typename TS> void
MultiVector<TV,TS>::insert(const std::vector<KeyValuePair> & data, SizeType smax){
	assert(smax>=0);
	assert(std::is_sorted(data.begin(),data.end(),
		[](KeyValuePair const &a, KeyValuePair const &b){return a.first<b.first;}));
	assert(!splits.empty());
	
	if(data.empty()) return;
	else {assert(smax>=data.back().first);}
				 
	splits.reserve(smax+1);
	values.reserve(values.size()+data.size());

	SizeType slast=data.front().first;
	assert(slast>=size());
	splits.insert(splits.end(), slast-size(), values.size());

	for(const auto [s,val] : data){
		if(slast!=s){
			splits.insert(splits.end(), s-slast, values.size());
			slast = s;
		}
		values.push_back(val);
	}
	
	if(smax==0) return;
	assert(smax>=slast);
	splits.insert(splits.end(),smax-slast, values.size());
}




#endif /* MultiVector_h */
