#include <math.h>

void reascl(float **a, int n, int n1, int n2)
{
	int i, j;
	float s;

	for (j=n1; j<=n2; j++) {
		s=0.0;
		for (i=1; i<=n; i++)
			if (fabs(a[i][j]) > fabs(s)) s=a[i][j];
		if (s != 0.0)
			for (i=1; i<=n; i++) a[i][j] /= s;
	}
}
