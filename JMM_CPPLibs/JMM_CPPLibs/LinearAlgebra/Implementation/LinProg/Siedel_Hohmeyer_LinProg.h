//
//  Siedel_Hohmeyer_LinProg.h
//  This header collects the files for
//  Michael Hohmeyer's implementation of R. Siedel's linear programming solver

// A few, very minor, adjustments had to be made to the source.

// Original source :
// https://github.com/cgrushko/seidel-lp-solver

// Academic paper :
// Seidel, R. (1991), "Small-dimensional linear programming and convex hulls made easy", Discrete & Computational Geometry 6 (1): 423–434, doi:10.1007/BF02574699

//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 15/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Siedel_Hohmeyer_LinProg_h
#define Siedel_Hohmeyer_LinProg_h

#define DOUBLE

#include "lp.h"
#include "localmath.h"
#include "tol.h"


#include "lp_base_case.c"
#include "linprog.c"
#include "vector_up.c"
#include "unit2.c"

#endif /* Siedel_Hohmeyer_LinProg_h */
