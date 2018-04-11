#include <math.h>

void besspq0(float x, float *p0, float *q0)
{
	if (x < 8.0) {
		float bessj0(float);
		void bessy01(float, float *, float *);
		float b,cosx,sinx,j0x,y0;
		b=sqrt(x)*1.25331413731550;
		bessy01(x,&y0,&j0x);
		j0x=bessj0(x);
		x -= 0.785398163397448;
		cosx=cos(x);
		sinx=sin(x);
		*p0 = b*(y0*sinx+j0x*cosx);
		*q0 = b*(y0*cosx-j0x*sinx);
	} else {
		int i;
		float x2,b0,b1,b2,y;
		static float ar1[11]={-0.10012e-15, 0.67481e-15, -0.506903e-14,
			0.4326596e-13, -0.43045789e-12, 0.516826239e-11,
			-0.7864091377e-10, 0.163064646352e-8, -0.5170594537606e-7,
			0.30751847875195e-5, -0.536522046813212e-3};
		static float ar2[10]={-0.60999e-15, 0.425523e-14,
			-0.3336328e-13, 0.30061451e-12, -0.320674742e-11,
			0.4220121905e-10, -0.72719159369e-9, 0.1797245724797e-7,
			-0.74144984110606e-6, 0.683851994261165e-4};
		y=8.0/x;
		x=2.0*y*y-1.0;
		x2=x+x;
		b1=b2=0.0;
		for (i=0; i<=10; i++) {
			b0=x2*b1-b2+ar1[i];
			b2=b1;
			b1=b0;
		}
		*p0 = x*b1-b2+0.99946034934752;
		b1=b2=0.0;
		for (i=0; i<=9; i++) {
			b0=x2*b1-b2+ar2[i];
			b2=b1;
			b1=b0;
		}
		*q0 = (x*b1-b2-0.015555854605337)*y;
	}
}
