#include <stdio.h>
void main ()
{
	void comdiv(float, float, float, float, float *, float *);
	float r,i;

	comdiv(-0.05,0.1,0.1,0.2,&r,&i);
	printf("(-.05+.1*i)/(.1+.2*i) = %-4.2f+%-4.2f*i",r,i);
}

