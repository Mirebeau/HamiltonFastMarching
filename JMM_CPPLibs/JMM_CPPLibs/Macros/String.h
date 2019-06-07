//
//  StringMacro.h
//  
//
//  Created by Jean-Marie Mirebeau on 21/09/2018.
//

#ifndef StringMacro_h
#define StringMacro_h

// Converts a preprocessor key into a string
#define STRING_NX(s) #s
#define STRING(s) STRING_NX(s)

#endif /* StringMacro_h */
