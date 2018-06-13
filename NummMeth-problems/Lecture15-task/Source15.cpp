#include <iostream>
#include <math.h>
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/tridag.h"
#include "../Source Code/code/utilities.h"
#include "../Lecture13-task/psplot.h"
using namespace std;

struct xi {
	Doub a, h;
	xi(Doub aa, Doub hh) :a(aa), h(hh) {}
	Doub operator() (Doub i)
	{
		return a + i * h;
	}
};


Doub unot(Doub x)
{
	return pow(x, 4.0);
}


int main() {
	// Initialize System
	Doub a = 0, b = 1.0;
	Doub dx = 0.4;
	Doub dt = dx;
	Doub D = 0.1;
	Doub t = 0.0;
	int n = (b - a) / dx - 1;
	int kkmax = 10;
	Doub tmax = 2.0;
	int m = tmax / dt;
	xi x(a, dx);
	VecDoub unow;
	VecDoub unext;
	VecDoub Va(n, 0.0); //Note that Jp is one shorter than J
	VecDoub Vb(n, 0.0);
	VecDoub Vc(n, 0.0); //Note that Jm is one shorter than J
	VecDoub F(n, 0.0);
	VecDoub uhalf(kkmax, 0.0);
	Doub s, olds = 0.0, oldolds = 0.0, rn = 0.0, r = 0.0, oldr = 0.0, oldoldr = 0.0;
	int kk = 0;
	while (kk<kkmax) {
		dx = dx / 2;
		//cout << "dx = " << dx << endl;
		dt = dx;
		n = (b - a) / dx - 1;
		m = tmax / dt;
		t = 0.0;
		xi xtemp(a, dx);
		VecDoub Temp(n, 0.0);
		x = xtemp;
		unow = Temp;
		unext = Temp;
		Va = Temp; //Note that Jp is one shorter than J
		Vb = Temp;
		Vc = Temp; //Note that Jm is one shorter than J
		F = Temp;

		//Init unow
		for (Doub i = 0; i<n; i++) {
			unow[i] = unot(x(i + 1));
		}

		//Time LOOP STARTS HERE
		for (int j = 0; j<m; j++) {
			// Define F vector
			int i = 0;
			F[i] = D / 2 * (unow[i + 1] - 2 * unow[i] + 0.0 + 0.0) / dx / dx + unow[i] / dt;
			i = n - 1;
			F[i] = D / 2 * (1.0 + 1.0 - 2 * unow[i] + unow[i - 1]) / dx / dx + unow[i] / dt;
			for (int i = 1; i<n - 1; i++) { F[i] = D / 2 * (unow[i + 1] - 2 * unow[i] + unow[i - 1]) / dx / dx + unow[i] / dt; }

			// Define Tridag Vectors
			for (int i = 0; i<n; i++)
			{
				Va[i] = -D / 2 / dx / dx;
				Vb[i] = D / dx / dx + 1.0 / dt;
				Vc[i] = -D / 2 / dx / dx;
			}
			// if you want, put them in a Matrix.. alternatively use Tridag
			//try {
				//tridag(Va, Vb, Vc, F, unext);
			//}
			//catch (NRerror &e) { cout << e.message << endl; }
			unow = unext;
			t += dt;

		}  // Time loop ends here

		   //util::print(unow,"unow");
		   //util::print(Va,"Va");
		   //util::print(Vb,"Vb");
		   //util::print(Vc,"Vc");
		   //util::print(F,"F");
		   //util::print(unext,"unext");
		cout << "x = " << x(n / 2 + 1) << setw(10);
		cout << "t = " << t << setw(10);
		uhalf[kk] = unow[n / 2];
		s = unow[n / 2];
		cout << setw(15) << kk + 1;
		cout << setw(15) << s;
		if (oldolds > 0.0) { cout << setw(15) << log2((oldolds - olds) / (olds - s)); }
		else { cout << setw(15) << " "; }
		if (olds > 0.0) { r = s + (s - olds) / 3; cout << setw(15) << r; } // We see that k1 tends towards 2,and with alpha=2 we get alpha^k-1 = 3
		else { cout << setw(15) << " "; }
		if (oldoldr > 0.0) { cout << setw(15) << log2((oldoldr - oldr) / (oldr - r)); }
		else { cout << setw(15) << " "; }
		if (oldr > 0.0) { rn = r + (r - oldr) / 15; cout << setw(15) << rn; } // We see that k2 tends towards 4,and with alpha=2 we get alpha^k-1 = 15
		else { cout << setw(15) << " "; }
		cout << endl;
		Doub alphak = abs((oldolds - olds) / (olds - s));

		if (kk > 3)
			if (abs(s - olds) / alphak < 0.0001)
				break;

		oldoldr = oldr;
		oldolds = olds;
		oldr = r;
		olds = s;
		kk++;
	} // kk loop ends here

	PSpage pg("myplot.ps");
	PSplot plot1(pg, 100., 500., 500., 800.);
	//	PSplot plot2(pg,100.,500.,100.,400.);
	plot1.setlimits(a, b, 0., 1.);
	plot1.frame();
	plot1.autoscales();
	plot1.xlabel("x");
	plot1.ylabel("y");
	//	plot1.lineplot(xlag,ylag);
	for (int i = 0; i< n + 2; i++)
	{
		plot1.pointsymbol(x(i), unot(x(i)), 108, 8);
	}
	for (int i = 0; i< unow.size(); i++)
	{
		plot1.pointsymbol(x(i + 1), unow[i], 105, 8);
	}
	pg.close();
	pg.display("\"ghostscript\" ");
}