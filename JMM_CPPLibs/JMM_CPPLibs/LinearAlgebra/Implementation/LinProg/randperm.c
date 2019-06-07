/* 
 * randperm.c
 *
 * Copyright (c) 1990 Michael E. Hohmeyer,
 *       hohmeyer@icemcfd.com
 * Permission is granted to modify and re-distribute this code in any manner
 * as long as this notice is preserved.  All standard disclaimers apply.
 *       
 */
int rand();
void randperm(int n,int perm[])
{
	int i, j, t;

	for(i=0; i<n; i++)
		perm[i] = i;
	for(i=0; i<n; i++) {
		j = rand()%(n-i)+i;
		t = perm[j];
		perm[j] = perm[i];
		perm[i] = t;
	}
}
