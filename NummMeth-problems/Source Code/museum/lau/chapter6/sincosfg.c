#include <math.h>

void sincosfg(float x, float *f, float *g)
{
	void sincosint(float, float *, float *);
	float chepolsum(int, float, float []);
	float absx,si,ci;

	absx=fabs(x);
	if (absx <= 4.0) {
		float cx,sx;
		sincosint(x,&si,&ci);
		cx=cos(x);
		sx=sin(x);
		si -= 1.570796326794897;
		*f = ci*sx-si*cx;
		*g = -ci*cx-si*sx;
	} else {
		float a[24];
		a[0] =  9.6578828035185e-1;  a[1]  = -4.3060837778597e-2;
		a[2] = -7.3143711748104e-3;  a[3]  =  1.4705235789868e-3;
		a[4] = -9.8657685732702e-5;  a[5]  = -2.2743202204655e-5;
		a[6] =  9.8240257322526e-6;  a[7]  = -1.8973430148713e-6;
		a[8] =  1.0063435941558e-7;  a[9]  =  8.0819364822241e-8;
		a[10]= -3.8976282875288e-8;  a[11] =  1.0335650325497e-8;
		a[12]= -1.4104344875897e-9;  a[13] = -2.5232078399683e-10;
		a[14]=  2.5699831325961e-10; a[15] = -1.0597889253948e-10;
		a[16]=  2.8970031570214e-11; a[17] = -4.1023142563083e-12;
		a[18]= -1.0437693730018e-12; a[19] =  1.0994184520547e-12;
		a[20]= -5.2214239401679e-13; a[21] =  1.7469920787829e-13;
		a[22]= -3.8470012979279e-14;
		*f = chepolsum(22,8.0/absx-1.0,a)/x;
		a[0] =  2.2801220638241e-1;  a[1]  = -2.6869727411097e-2;
		a[2] = -3.5107157280958e-3;  a[3]  =  1.2398008635186e-3;
		a[4] = -1.5672945116862e-4;  a[5]  = -1.0664141798094e-5;
		a[6] =  1.1170629343574e-5;  a[7]  = -3.1754011655614e-6;
		a[8] =  4.4317473520398e-7;  a[9]  =  5.5108696874463e-8;
		a[10]= -5.9243078711743e-8;  a[11] =  2.2102573381555e-8;
		a[12]= -5.0256827540623e-9;  a[13] =  3.1519168259424e-10;
		a[14]=  3.6306990848979e-10; a[15] = -2.2974764234591e-10;
		a[16]=  8.5530309424048e-11; a[17] = -2.1183067724443e-11;
		a[18]=  1.7133662645092e-12; a[19] =  1.7238877517248e-12;
		a[20]= -1.2930281366811e-12; a[21] =  5.7472339223731e-13;
		a[22]= -1.8415468268314e-13; a[23] =  3.5937256571434e-14;
		*g = 4.0*chepolsum(23,8.0/absx-1.0,a)/absx/absx;
	}
}
