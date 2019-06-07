/* 
 * lp.h
 *
 * Copyright (c) 1990 Michael E. Hohmeyer, 
 *       hohmeyer@icemcfd.com
 * Permission is granted to modify and re-distribute this code in any manner
 * as long as this notice is preserved.  All standard disclaimers apply.
 *
 * version history
 *     1/22/1991 - original version
 *     1/5/1995  - fix bug in projection of degenerate objective
 *                 function to lower dimension
 */
/* status from lp_intersect or linprog */ 
#define INFEASIBLE 0
#define MINIMUM 1
#define UNBOUNDED 2
#define AMBIGUOUS 3

/* status from plane_down */
#define REDUNDANT 0
#define PROPER 1

#include "tol.h"

#ifdef DOUBLE
#define linprog(v,istart, n,num,den,dim,opt,work,next,prev,max_size)  \
dlinprog(v,istart, n,num,den,dim,opt,work,next,prev,max_size)
#else
#define linprog(v,istart, n,num,den,dim,opt,work,next,prev,max_size)  \
slinprog(v,istart, n,num,den,dim,opt,work,next,prev,max_size)
#endif

/*
#ifdef __cplusplus
extern "C" {
#endif
 */
    void randperm(int i, int *p);
    void randomize(int n, int *perm);
    int linprog(FLOAT *v, int istart,int n, FLOAT *num, FLOAT *den,
        int dim, FLOAT *opt, FLOAT *work, int *next, int *prev, int max_size);
    int lp_base_case(FLOAT halves[][2], int m, FLOAT n_vec[2],
        FLOAT d_vec[2], FLOAT opt[2], int *next, int *prev);
    int wedge(FLOAT halves[][2], int m, int *next, int *prev,
        FLOAT cw_vec[2], FLOAT ccw_vec[2], int *degen, int max_halves);
    void plane_down(FLOAT *elim_eqn, int ivar, int idim,
        FLOAT *old_plane, FLOAT *new_plane);
    void findimax(FLOAT *pl,int idim,int *imax);
    void vector_up(FLOAT *equation,int ivar,int idim,
        FLOAT *low_vector,FLOAT *vec);
/*
#ifdef __cplusplus
}
#endif
 */
