#ifndef TOL_H
#define TOL_H

#ifdef DOUBLE

/* this is the worst error we'll get doing a few additions
 * of numbers on the order of unity */
#define EPS		2.0e-16
#define FLOAT	double

#else

#define EPS 		1.0e-7
#define FLOAT	float

#endif

#endif
