#ifndef AD2_hpp
#define AD2_hpp

template<typename TAD, typename TComponent>
struct AD2Operators {
	typedef TAD AD;
	typedef TComponent Comp;
	friend AD operator + (AD a, AD const & b) { a += b; return a; }
	friend AD operator + (AD a, Comp const & b) { a += b; return a; }
	friend AD operator + (Comp const & a, AD b) { b += a; return b; }
	
	friend AD operator - (AD a, AD const & b) { a -= b; return a; }
	friend AD operator - (AD a, Comp const & b) { a -= b; return a; }
	friend AD operator - (Comp const & a, AD b) { b -= a; return -b; }
	
	friend AD operator * (AD a, AD const & b) { a *= b; return a; }
	friend AD operator * (AD a, Comp const & b) { a *= b; return a; }
	friend AD operator * (Comp const & a, AD b) { b *= a; return b; }
};

#endif
