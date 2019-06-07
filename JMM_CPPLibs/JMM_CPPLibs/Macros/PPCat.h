//
//  PPCatMacro.h
//  
//
//  Created by Jean-Marie Mirebeau on 21/09/2018.
//

#ifndef PPCatMacro_h
#define PPCatMacro_h

// Concatenates two preprocessor keys
#define PPCAT_NX(A, B) A ## B
#define PPCAT(A, B) PPCAT_NX(A, B)

#endif /* PPCatMacro_h */
