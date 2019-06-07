//
//  TemplateLog2.h
//  FileHFM
//
//  Created by Jean-Marie Mirebeau on 16/11/2018.
//

#ifndef TemplateLog2_h
#define TemplateLog2_h

// Silly replacement for a constexpr log. TODO : hide somewhere.

template<size_t i> struct FloorLog2 {static const size_t value = 1+FloorLog2<i/2>::value;};
template<> struct FloorLog2<1> {static const size_t value = 0;};
template<size_t i> class CeilLog2 {
	static const size_t floor = FloorLog2<i>::value;
public:
	static const size_t value = (i== (1<<floor) ? floor : (floor+1));
};

#endif /* TemplateLog2_h */
