#include <math.h>

void besspq1(float x, float *p1, float *q1)
{
	if (x < 8.0) {
		float bessj1(float);
		void bessy01(float, float *, float *);
		float b,cosx,sinx,j1x,y1;
		b=sqrt(x)*1.25331413731550;
		bessy01(x,&j1x,&y1);
		j1x=bessj1(x);
		x -= 0.785398163397448;
		cosx=cos(x);
		sinx=sin(x);
		*p1 = b*(j1x*sinx-y1*cosx);
		*q1 = b*(j1x*cosx+y1*sinx);
	} else {
		int i;
		float x2,b0,b1,b2,y;
		static float ar1[11]={0.10668e-15, -0.72212e-15, 0.545267e-14,
			-0.4684224e-13, 0.46991955e-12, -0.570486364e-11,
			0.881689866e-10, -0.187189074911e-8, 0.6177633960644e-7,
			-0.39872843004889e-5, 0.89898983308594e-3};
		static float ar2[11]={-0.10269e-15, 0.65083e-15, -0.456125e-14,
			0.3596777e-13, -0.32643157e-12, 0.351521879e-11,
			-0.4686363688e-10, 0.82291933277e-9, -0.2095978138408e-7,
			0.91386152579555e-6, -0.96277235491571e-4};
		y=8.0/x;
		x=2.0*y*y-1.0;
		x2=x+x;
		b1=b2=0.0;
		for (i=0; i<=10; i++) {
			b0=x2*b1-b2+ar1[i];
			b2=b1;
			b1=b0;
		}
		*p1 = x*b1-b2+1.0009030408600137;
		b1=b2=0.0;
		for (i=0; i<=10; i++) {
			b0=x2*b1-b2+ar2[i];
			b2=b1;
			b1=b0;
		}
		*q1 = (x*b1-b2+0.46777787069535e-1)*y;
	}
}
