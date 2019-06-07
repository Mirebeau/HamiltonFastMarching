/* 
 * vector_up.c
 *
 * Copyright (c) 1990 Michael E. Hohmeyer,
 *       hohmeyer@icemcfd.com
 * Permission is granted to modify and re-distribute this code in any manner
 * as long as this notice is preserved.  All standard disclaimers apply.
 *       
 */
#include "tol.h"
#include "lp.h"


void vector_up(FLOAT equation[],int ivar,int idim,
		FLOAT low_vector[],FLOAT vector[])
{
	int i;

	vector[ivar] = 0.0;
	for(i=0; i<ivar; i++) {
		vector[i] = low_vector[i];
		vector[ivar] -= equation[i]*low_vector[i];
	}
	for(i=ivar+1; i<=idim; i++) {
		vector[i] = low_vector[i-1];
		vector[ivar] -= equation[i]*low_vector[i-1];
	}
	vector[ivar] /= equation[ivar];
}
