#include <math.h>

#ifdef DOUBLE
int
d_unit2(double a[],double b[],double eps)
{
	double size;
	size = sqrt(a[0]*a[0] + a[1]*a[1]);
	if(size < 2*eps) return(1);
	b[0] = a[0]/size;
	b[1] = a[1]/size;
	return(0);
}
#else
int
s_unit2(float a[],float b[],float eps)
{
	float size;
	size = sqrt(a[0]*a[0] + a[1]*a[1]);
	if(size < 2*eps) return(1);
	b[0] = a[0]/size;
	b[1] = a[1]/size;
	return(0);
}
#endif
