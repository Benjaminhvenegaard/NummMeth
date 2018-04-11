#include <math.h>

void hsh2row3(int l, int ua, int ub, int ux, int j, float a1, float a2,
					float **a, float **b, float **x)
{
	float *allocate_real_vector(int, int);
	void free_real_vector(float *, int);
	void hshvectam(int, int, int, int, float, float [], float **);
	float *v,d1,d2,s1,s2,r,d,c;

	if (a2 != 0.0) {
		v=allocate_real_vector(j,j+1);
		d1=fabs(a1);
		d2=fabs(a2);
		s1 = (a1 >= 0.0) ? 1.0 : -1.0;
		s2 = (a2 >= 0.0) ? 1.0 : -1.0;
		if (d2 <= d1) {
			r=d2/d1;
			d=sqrt(1.0+r*r);
			c = -1.0-1.0/d;
			v[j]=s1*s2*r/(1.0+d);
		} else {
			r=d1/d2;
			d=sqrt(1.0+r*r);
			c = -1.0-r/d;
			v[j]=s1*s2/(r+d);
		}
		v[j+1]=1.0;
		hshvectam(l,ua,j,j+1,c,v,a);
		hshvectam(l,ub,j,j+1,c,v,b);
		hshvectam(1,ux,j,j+1,c,v,x);
		free_real_vector(v,j);
	}
}
