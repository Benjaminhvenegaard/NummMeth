#include <stdio.h>
void main ()
{
	void bessjaplusn(float, float, int, float []);
	int n;
	float a,x,ja[3];

	x=2.0;
	a=0.78;
	n=2;
	bessjaplusn(a,x,n,ja);
	printf("BESSJAPLUSN delivers:\n"
			" X = %2.0f    A = %4.2f    N = %d\n"
			" %e   %e   %e\n",x,a,n,ja[0],ja[1],ja[2]);
}

