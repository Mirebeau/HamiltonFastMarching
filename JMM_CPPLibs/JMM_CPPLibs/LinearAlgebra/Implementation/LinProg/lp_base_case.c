/* 
 * lp_base_case.c
 *
 * Copyright (c) 1990 Michael E. Hohmeyer,
 *       hohmeyer@icemcfd.com
 * Permission is granted to modify and re-distribute this code in any manner
 * as long as this notice is preserved.  All standard disclaimers apply.
 *
 */
#include <math.h>
#include "localmath.h"
#define not_zero(a) ((a) >= 2*EPS || (a) <= -2*EPS)

#include "lp.h"

int lp_no_constraints(int d,FLOAT n_vec[],FLOAT d_vec[],FLOAT opt[]);

void lp_min_lin_rat(int degen, FLOAT cw_vec[2], 
		FLOAT ccw_vec[2], FLOAT n_vec[2], FLOAT d_vec[2], FLOAT opt[2]);

int move_to_front(int i,int next[],int prev[],int max_halves);


/*
int
unit2(FLOAT a[],FLOAT b[],FLOAT eps)
{
	FLOAT size;

	size = fsqrt(a[0]*a[0] + a[1]*a[1]);
	if(size < 2*eps) return(1);
	b[0] = a[0]/size;
	b[1] = a[1]/size;
	return(0);
}
*/
/*
 * return the minimum on the projective line
 *
 */
int lp_base_case(FLOAT halves[][2], 	/* halves --- half lines */
	int m, 				/* m      --- terminal marker */
	FLOAT n_vec[2],			/* n_vec  --- numerator funciton */
	FLOAT d_vec[2],			/* d_vec  --- denominator function */
	FLOAT opt[2],			/* opt    --- optimum  */
	int next[],
	int prev[],
    int max_halves)			/* next, prev  ---
					double linked list of indices */
{
	FLOAT cw_vec[2], ccw_vec[2];
#ifdef CHECK
	FLOAT d_cw;
	int i;
#endif
	int status, degen;
	FLOAT ab;

/* find the feasible region of the line */
	status = wedge(halves,m,next,prev,cw_vec,ccw_vec,&degen,max_halves);

	if(status==INFEASIBLE) return(status);
/* no non-trivial constraints one the plane: return the unconstrained
** optimum */
	if(status==UNBOUNDED) {
		return(lp_no_constraints(1,n_vec,d_vec,opt));
	}
        ab = ABS(cross2(n_vec,d_vec));
	if(ab < 2*EPS*EPS) {
		if(dot2(n_vec,n_vec) < 2*EPS*EPS ||
			 dot2(d_vec,d_vec) > 2*EPS*EPS) {
/* numerator is zero or numerator and denominator are linearly dependent */
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
			status = AMBIGUOUS;
		} else {
/* numerator is non-zero and denominator is zero 
** minimize linear functional on circle */
			if(!degen && cross2(cw_vec,n_vec) <= 0.0 &&
			cross2(n_vec,ccw_vec) <= 0.0 ) {
/* optimum is in interior of feasible region */
				opt[0] = -n_vec[0];
				opt[1] = -n_vec[1];
			} else if(dot2(n_vec,cw_vec) > dot2(n_vec,ccw_vec) ) {
/* optimum is at counter-clockwise boundary */
				opt[0] = ccw_vec[0];
				opt[1] = ccw_vec[1];
			} else {
/* optimum is at clockwise boundary */
				opt[0] = cw_vec[0];
				opt[1] = cw_vec[1];
			}
			status = MINIMUM;
		}
	} else {
/* niether numerator nor denominator is zero */
		lp_min_lin_rat(degen,cw_vec,ccw_vec,n_vec,d_vec,opt);
		status = MINIMUM;
	}
#ifdef CHECK
	for(i=0; i!=m; i=next[i]) {
		d_cw = dot2(opt,halves[i]);
		if(d_cw < -2*EPS) {
			printf("error at base level\n");
			exit(1);
		}
	}
#endif
	return(status);
}
void lp_min_lin_rat(int degen,
		FLOAT cw_vec[2],
		FLOAT ccw_vec[2],
		FLOAT n_vec[2],
		FLOAT d_vec[2],
		FLOAT opt[2])
{
	FLOAT d_cw, d_ccw, n_cw, n_ccw;

/* linear rational function case */
	d_cw = dot2(cw_vec,d_vec);
	d_ccw = dot2(ccw_vec,d_vec);
	n_cw = dot2(cw_vec,n_vec);
	n_ccw = dot2(ccw_vec,n_vec);
	if(degen) {
/* if degenerate simply compare values */
		if(n_cw/d_cw < n_ccw/d_ccw) {
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
/* check that the clock-wise and counter clockwise bounds are not near a poles */
	} else if(not_zero(d_cw) && not_zero(d_ccw)) {
/* the valid region does not contain a poles */
		if(d_cw*d_ccw > 0.0) {
/* find which end has the minimum value */
			if(n_cw/d_cw < n_ccw/d_ccw) {
				opt[0] = cw_vec[0];
				opt[1] = cw_vec[1];
			} else {
				opt[0] = ccw_vec[0];
				opt[1] = ccw_vec[1];
			}
		} else {
/* the valid region does contain a poles */
			if(d_cw > 0.0) {
				opt[0] = -d_vec[1];
				opt[1] = d_vec[0];
			} else {
				opt[0] = d_vec[1];
				opt[1] = -d_vec[0];
			}
		}
	} else if(not_zero(d_cw)) {
/* the counter clockwise bound is near a pole */
		if(n_ccw*d_cw > 0.0) {
/* counter clockwise bound is a positive pole */
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
/* counter clockwise bound is a negative pole */
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
	} else if(not_zero(d_ccw)) {
/* the clockwise bound is near a pole */
		if(n_cw*d_ccw > 2*EPS) {
/* clockwise bound is at a positive pole */
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		} else {
/* clockwise bound is at a negative pole */
			 opt[0] = cw_vec[0];
			 opt[1] = cw_vec[1];
		} 
	} else {
/* both bounds are near poles */
		if(cross2(d_vec,n_vec) > 0.0) {
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
	}
}
int wedge(FLOAT halves[][2],
	int m,
	int next[],
	int prev[],
	FLOAT cw_vec[],
	FLOAT ccw_vec[],
	int *degen,
    int max_halves)
{
	int i;
	FLOAT d_cw, d_ccw;
	int offensive;

	*degen = 0;
	for(i=0;i!=m;i = next[i]) {
		if(!unit2(halves[i],ccw_vec,EPS)) {
/* clock-wise */
			cw_vec[0] = ccw_vec[1];
			cw_vec[1] = -ccw_vec[0];
/* counter-clockwise */
			ccw_vec[0] = -cw_vec[0];
			ccw_vec[1] = -cw_vec[1];
			break;
		}
	}
	if(i==m) return(UNBOUNDED);
	i = 0;
	while(i!=m) {
		offensive = 0;
		d_cw = dot2(cw_vec,halves[i]);
		d_ccw = dot2(ccw_vec,halves[i]);
		if(d_ccw >= 2*EPS) {
			if(d_cw <= -2*EPS) {
				cw_vec[0] = halves[i][1];
				cw_vec[1] = -halves[i][0];
				(void)unit2(cw_vec,cw_vec,EPS);
				offensive = 1;
			}
		} else if(d_cw >= 2*EPS) {
			if(d_ccw <= -2*EPS) {
				ccw_vec[0] = -halves[i][1];
				ccw_vec[1] = halves[i][0];
				(void)unit2(ccw_vec,ccw_vec,EPS);
				offensive = 1;
			}
		} else if(d_ccw <= -2*EPS && d_cw <= -2*EPS) {
			return(INFEASIBLE);
		} else if((d_cw <= -2*EPS) ||
			(d_ccw <= -2*EPS) ||
			(cross2(cw_vec,halves[i]) < 0.0)) {
/* degenerate */
			if(d_cw <= -2*EPS) {
				(void)unit2(ccw_vec,cw_vec,EPS);
			} else if(d_ccw <= -2*EPS) { 
				(void)unit2(cw_vec,ccw_vec,EPS);
			}
			*degen = 1;
			offensive = 1;
		}
/* place this offensive plane in second place */
		if(offensive) i = move_to_front(i,next,prev,max_halves);
		i = next[i];
		if(*degen) break;
	}
	if(*degen) {
		while(i!=m) {
			d_cw = dot2(cw_vec,halves[i]);
			d_ccw = dot2(ccw_vec,halves[i]);
			if(d_cw < -2*EPS) {
				if(d_ccw < -2*EPS) {
					return(INFEASIBLE);
				} else {
					cw_vec[0] = ccw_vec[0];
					cw_vec[1] = ccw_vec[1];
				}
			} else if(d_ccw < -2*EPS) {
				ccw_vec[0] = cw_vec[0];
				ccw_vec[1] = cw_vec[1];
			}
			i = next[i];
		}
	}
	return(MINIMUM);
}
