#ifndef LOCALMATH
#define LOCALMATH

#include "tol.h"

typedef FLOAT twoD[2];
#define ABS(a) ((a)>0.0 ? (a) : -(a))

/* line dot product line vs. point */
#define ldot(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2])

/*
 * line cross product a line from two points
 *  (az, ay, 1) x (bx, bx, 1) 
 */
#define lcross(a, b, c) {       (c)[0] = (a)[1] - (b)[1]; \
				(c)[1] = (b)[0] - (a)[0]; \
				(c)[2] = (a)[0]*(b)[1] - (a)[1]*(b)[0]; }

#define vsub2(a, b, c) 	{ (c)[0] = (a)[0] - (b)[0]; \
                        (c)[1] = (a)[1] - (b)[1]; }
#define vadd2(a, b, c)	{ (c)[0] = (a)[0] + (b)[0]; \
			(c)[1] = (a)[1] + (b)[1]; }
#define vcopy2(a, b) 	{ (b)[0] = (a)[0]; \
                        (b)[1] = (a)[1]; }

#define vmult2(s,a,b) 	{ (b)[0] = (s)*(a)[0]; \
			(b)[1] = (s)*(a)[1]; }
#define dot2(a,b) 	((a)[0]*(b)[0] + (a)[1]*(b)[1])
#define cross2(a,b) 	((a)[0]*(b)[1] - (a)[1]*(b)[0])
#define aprxlen2(a)	(ABS((a)[0]) + ABS((a)[1]))


#define midpoint(a, b, c)	{ (c)[0] = ((a)[0]  + (b)[0])*.5; \
				 (c)[1] = ((a)[1]  + (b)[1])*.5; \
				 (c)[2] = ((a)[2]  + (b)[2])*.5; }

#define dot(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
//#define min(a,b) ((a)<(b) ?(a):(b)) // Commented by Mirebeau. Shadows std::max, std::min
//#define max(a,b) ((a)>(b) ?(a):(b))
#define square(a) ((a)*(a))
typedef FLOAT point[3];
typedef FLOAT vector[3];
typedef FLOAT triple[3];
typedef FLOAT plane[4];

#define vadd(a, b, c)	{ (c)[0] = (a)[0] + (b)[0]; \
			(c)[1] = (a)[1] + (b)[1]; \
			(c)[2] = (a)[2] + (b)[2]; }
#define vsub(a, b, c)	{ (c)[0] = (a)[0] - (b)[0]; \
			(c)[1] = (a)[1] - (b)[1]; \
			(c)[2] = (a)[2] - (b)[2]; }
#define vmult(s, a, b) 	{ (b)[0] = (s)*(a)[0]; \
			(b)[1] = (s)*(a)[1]; \
			(b)[2] =(s)*(a)[2]; }
#define size2(a) dot((a),(a))
#define vzero2(a)	{ (a)[0] = (a)[1] = 0.0; }


#define vlinc(s, a, t, b, c)  { (c)[0] = (s)*(a)[0] + (t)*(b)[0]; \
				(c)[1] = (s)*(a)[1] + (t)*(b)[1]; \
				(c)[2] = (s)*(a)[2] + (t)*(b)[2]; }

/* multiply and add */
#define vmadd(s, a, b, c)  { (c)[0] = (s)*(a)[0] + (b)[0]; \
				(c)[1] = (s)*(a)[1] + (b)[1]; \
				(c)[2] = (s)*(a)[2] + (b)[2]; }

#define vzero(a)	{ (a)[0] = (a)[1] = (a)[2] = 0.0; }

#define vcopy(a, b)	{ (b)[0] = (a)[0]; \
			  (b)[1] = (a)[1]; \
			  (b)[2] = (a)[2]; }
#define vneg(a,b) 	{ (b)[0] = -(a)[0]; \
			  (b)[1] = -(a)[1]; \
			  (b)[2] = -(a)[2]; }


#ifdef DOUBLE

#define unit2(a,b,eps)		d_unit2(a,b,eps)
#define triple_cross(a,b,c)	d_triple_cross(a,b,c)
#define cross(a,b,c)		d_cross(a,b,c)
#define unit(a,b)		d_unit(a,b)
#define solve3(a,x,b)		d_solve3(a,x,b)
#define gen_perp(a,b)		d_gen_perp(a,b)
#define distance(a,b)		d_distance(a,b)
#define vswap(a,b)		d_vswap(a,b)
#define Fsort(a,b)		d_sort(a,b)
#define Fmerge_sort(a,b,c)	d_merge_sort(a,b,c)

#else

#define unit2(a,b,eps)		s_unit2(a,b,eps)
#define triple_cross(a,b,c)	s_triple_cross(a,b,c)
#define cross(a,b,c)		s_cross(a,b,c)
#define unit(a,b)		s_unit(a,b)
#define solve3(a,x,b)		s_solve3(a,x,b)
#define gen_perp(a,b)		s_gen_perp(a,b)
#define distance(a,b)		s_distance(a,b)
#define vswap(a,b)		s_vswap(a,b)
#define Fsort(a,b)		s_sort(a,b)
#define Fmerge_sort(a,b,c)	s_merge_sort(a,b,c)

#endif


#ifdef __cplusplus
extern "C" {
#endif

void i_sort(int size, int *list);

#ifdef DOUBLE

int d_unit2(double *a, double *b, double eps);
double d_triple_cross(double *a, double *b, double *c);
void d_cross(double *a, double *b, double *c);
int d_unit(double *a, double *b); 
int d_solve3(double a[3][3], double *b, double *x);
void d_gen_perp(double *v, double *perp);
double d_distance(double *a, double *b);
void d_vswap(double *a, double *b);
void d_sort(int i, double *a);
void d_merge_sort(int size, float *list, float *temp);

#else 

int s_unit2(float *a, float *b, float eps);
float s_triple_cross(float *a, float *b, float *c);
void s_cross(float *a, float *b, float *c);
int s_unit(float *a, float *b); 
int s_solve3(float a[3][3], float *b, float *x);
void s_gen_perp(float *v, float *perp);
float s_distance(float *a, float *b);
void s_vswap(float *a, float *b);
void s_sort(int i, float *a);
void s_merge_sort(int size, float *list, float *temp);
int invert3(float from[3][3], float to[3][3]);

#endif


#ifdef __cplusplus
}
#endif

#define SWAP(type,a,b) {type temp = (a); (a) = (b); (b) = temp; }
#define ORDER(type, a, b)       if((a) > (b)) SWAP(type, a, b);

#endif
