// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef RedeclareTypesMacro_h
#define RedeclareTypesMacro_h

/* Usage : 
 #define FromA(x) A:: x

 struct A {
    typedef int Integer;
    typedef double Scalar;
    cont int d=2;
 };
 
 struct B {
    Redeclare2Types(FromA,Integer,Scalar)
    Redeclare1Constant(FromA,d)
 };
 
 // For constants, in combination with templates, the linker may unfortunately require one additional line such as
 constexpr typeof(A::d) B::d
 */

#include "../Macros/Expand.h"
#define insertdd(a,b) a :: b

#define Redeclare1Type(From,__a) typedef typename insertdd(From,__a) __a;
#define Redeclare2Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare1Type(From,__VA_ARGS__))
#define Redeclare3Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare2Types(From,__VA_ARGS__))
#define Redeclare4Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare3Types(From,__VA_ARGS__))
#define Redeclare5Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare4Types(From,__VA_ARGS__))
#define Redeclare6Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare5Types(From,__VA_ARGS__))
#define Redeclare7Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare6Types(From,__VA_ARGS__))
#define Redeclare8Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare7Types(From,__VA_ARGS__))
#define Redeclare9Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare8Types(From,__VA_ARGS__))
#define Redeclare10Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare9Types(From,__VA_ARGS__))
#define Redeclare11Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare10Types(From,__VA_ARGS__))
#define Redeclare12Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare11Types(From,__VA_ARGS__))
#define Redeclare13Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare12Types(From,__VA_ARGS__))
#define Redeclare14Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare13Types(From,__VA_ARGS__))
#define Redeclare15Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare14Types(From,__VA_ARGS__))
#define Redeclare16Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare15Types(From,__VA_ARGS__))
#define Redeclare17Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare16Types(From,__VA_ARGS__))
#define Redeclare18Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare17Types(From,__VA_ARGS__))
#define Redeclare19Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare18Types(From,__VA_ARGS__))
#define Redeclare20Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare19Types(From,__VA_ARGS__))
#define Redeclare21Types(From,__a,...) Redeclare1Type(From,__a) EXPAND(Redeclare20Types(From,__VA_ARGS__))


#define Redeclare1Constant(From,__a) static constexpr decltype(insertdd(From,__a)) __a = insertdd(From,__a);
#define Redeclare2Constants(From,__a,...) Redeclare1Constant(From,__a) EXPAND(Redeclare1Constant(From,__VA_ARGS__))
#define Redeclare3Constants(From,__a,...) Redeclare1Constant(From,__a) EXPAND(Redeclare2Constants(From,__VA_ARGS__))
#define Redeclare4Constants(From,__a,...) Redeclare1Constant(From,__a) EXPAND(Redeclare3Constants(From,__VA_ARGS__))
#define Redeclare5Constants(From,__a,...) Redeclare1Constant(From,__a) EXPAND(Redeclare4Constants(From,__VA_ARGS__))


#endif /* RedeclareTypesMacro_h */
