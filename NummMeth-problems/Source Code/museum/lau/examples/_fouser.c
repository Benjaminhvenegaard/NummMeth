#include <math.h>
#include <stdio.h>
void main ()
{
	float fouser(int, float, float []);
	float pi, a[2] = {0.5, 1.0};

	pi=atan(1.0)*4.0;
	printf("FOUSER delivers:  %-7.2f%-7.2f%-7.2f",
			fouser(1,0.0,a),fouser(1,pi/2.0,a),fouser(1,pi,a));
}

