#include <stdio.h>
void main ()
{
	float airyzeros(int, int, float [], float []);
	float b,zbi[4],vbid[4];

	b=airyzeros(3,2,zbi,vbid);
	printf("AIRYZEROS delivers:\n"
			" The third zero of BI(X) is  %e\n"
			" The value of (D/DX)BI(X) in this point is  %e\n",
			b,vbid[3]);
}

