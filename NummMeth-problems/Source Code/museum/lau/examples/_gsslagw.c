#include <math.h>
#include <stdio.h>
void main ()
{
	void gsslagwghts(int, float, float [], float []);
	int n;
	float ind,x[11],w[11];

	gsslagwghts(10,0.0,x,w);
	ind=0.0;
	for (n=10; n>=1; n--) ind += w[n]*sin(x[n]);
	printf("Delivers:  %13.6e\n",ind-0.5);
}

