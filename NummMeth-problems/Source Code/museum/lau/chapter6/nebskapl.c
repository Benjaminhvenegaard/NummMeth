#include <math.h>

void nonexpbesskaplusn(float a, float x, int nmax, float kan[])
{
	void nonexpbesska01(float, float, float *, float *);
	int n;
	float k1;

	nonexpbesska01(a,x,&kan[0],&k1);
	a -= 1.0;
	x=2.0/x;
	if (nmax > 0) kan[1]=k1;
	for (n=2; n<=nmax; n++) kan[n]=kan[n-2]+(a+n)*x*kan[n-1];
}
