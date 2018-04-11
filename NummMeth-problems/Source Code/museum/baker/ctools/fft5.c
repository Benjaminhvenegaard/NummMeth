/*  fast Fourier transform (FFT)


from Handbook of C tools for Scientists and Engineers by L. Baker

DEPENDENCIES:

header file complex.h required

*/
#include <stdio.h>
#include "complex.h"
#define twopi 6.283185307

int iterp;/* global to return count used*/

#define max(a,b) (((a)>(b))? (a): (b))
#define abs(x) ((x)?  (x):-(x))
#define DOFOR(i,to) for(i=0;i<to;i++)

main(argc,argv) int argc;char **argv;
{int i,ii,nh,n=16;
double invn,exp(),dt=.25,omega,realpt,imagpt;
struct complex w[32],wi[32], data[32];
nh=n>>1;
DOFOR(i,n)
	{
	CMPLX(data[i], exp(-dt*(i)),0.);
	};
/* caveat*/ data[0].x=.5;/*not 1. see text*/
printf(" before transform:\n");
DOFOR(i,n){printc(&(data[i])) ;printf("\n");}
invn=1./n;
fftinit(w,wi,n);
printf(" w factors:\n");
DOFOR(i,n){printc(&(w[i])) ;printc(&(wi[i])) ;printf("\n");}
printf(" transformed:\n");
fft(data,w,n);
DOFOR(i,n){printc(&(data[i])) ;printf("\n");}
printf(" scaled by dt and compared to analytic answer:\n");
DOFOR(i,n)
	{
	ii= (i-nh)<0 ? i : i-n ;
	omega= twopi*ii/(n*dt);
	realpt= 1./(1.+omega*omega);
	imagpt= -realpt*omega;
	printf(" %f %f  %f %f\n", dt*data[i].x,dt*data[i].y,realpt,imagpt);
	}
/* inverse fft*/
fft(data,wi,n);
printf(" transformed back and scaled by 1/N:\n");
DOFOR(i,n)
	{
	CTREAL(data[i],data[i],(invn));
	printc(&(data[i]));printf("\n");
	}
exit(0);
}



int bitr(k,logn) int k,logn;
{
int ans,j,i;
ans=0;
j=k;
DOFOR(i,logn)
	{
	ans=(ans<<1)+(j&1);
	j=j>>1;
	}
return(ans);
}


int log2(n) int n;
{
int i;
i=-1;/* will return -1 if n<=0 */
while(1)
	{
	if(n==0)break;
	n=n>>1;
	i++;
	}
return(i);
}

printc(x) struct complex *x;
{
printf("%f %f",x->x,x->y);
return;
}

fftinit(w,wi,n) int n; struct complex w[],wi[];
{
int i;
double realpt,imagpt,cos(),sin();
double factr,angle;
factr=twopi/n;
DOFOR(i,n)
	{angle=i*factr;
	 realpt=cos(angle);imagpt=sin(angle);
	CMPLX(w[i],realpt,imagpt);
	CMPLX(wi[i],realpt,(-imagpt));
	}
return;
}


fft(x,w,n) int n; struct complex w[],x[];
{
int n1,logn,i,j,k,l,logl,p;
struct complex s,t;
logn=log2(n);
n1=n>>1;
j=logn-1;
/* transform*/
k=0;
DOFOR(logl,logn)
	{
	do{
	DOFOR(i,n1)
		{
		p=bitr((k>>j),logn);
		l=k+n1;
			CONJG(s,w[p]);
			CMULT(t,s,x[l]);
			CSUB(x[l],x[k],t);
			CADD(x[k],t,x[k]);
			k++;
		};/* dofor i*/
	k+=n1;
	}while (k<n);
k=0;
j--;
n1=n1>>1;
}


/*reorder*/
for(i=1;i<n;i++)
	{
	k=bitr(i,logn);
	if (i>k)
		{
		/*exchange i,k elements*/
		CLET(s,x[i]);
		CLET(x[i],x[k]);
		CLET(x[k],s);		
		}
	
	};

return;
}

