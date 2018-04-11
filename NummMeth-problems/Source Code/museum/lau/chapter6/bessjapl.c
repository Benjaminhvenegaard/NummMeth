#include <math.h>

void bessjaplusn(float a, float x, int n, float ja[])
{
	if (x == 0.0) {
		ja[0] = (a == 0.0) ? 1.0 : 0.0;
		for (; n>=1; n--) ja[n]=0.0;
	} else if (a == 0.0) {
		void bessj(float, int, float []);
		bessj(x,n,ja);
	} else if (a == 0.5) {
		void spherbessj(float, int, float []);
		float s;
		s=sqrt(x)*0.797884560802865;
		spherbessj(x,n,ja);
		for (; n>=0; n--) ja[n] *= s;
	} else {
		float gamma(float);
		int start(float, int, int);
		int k,m,nu;
		float a2,x2,r,s,l,labda;
		l=1.0;
		nu=start(x,n,0);
		for (m=1; m<=nu; m++) l=l*(m+a)/(m+1);
		r=s=0.0;
		x2=2.0/x;
		k = -1;
		a2=a+a;
		for (m=nu+nu; m>=1; m--) {
			r=1.0/(x2*(a+m)-r);
			if (k == 1)
				labda=0.0;
			else {
				l=l*(m+2)/(m+a2);
				labda=l*(m+a);
			}
			s=r*(labda+s);
			k = -k;
			if (m <= n) ja[m]=r;
		}
		ja[0]=r=1.0/gamma(1.0+a)/(1.0+s)/pow(x2,a);
		for (m=1; m<=n; m++) r = ja[m] *= r;
	}
}
