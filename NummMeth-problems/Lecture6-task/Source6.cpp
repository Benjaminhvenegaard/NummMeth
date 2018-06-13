#include <iostream>
#include <math.h>


#include "../Source Code/code/nr3.h"
#include "../Source Code/code/ludcmp.h"
#include "../Source Code/code/qrdcmp.h"
//#include <roots_multidim.h>

using namespace std;

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p, VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
	const Doub ALF = 1.0e-4, TOLX = numeric_limits<Doub>::epsilon();//te
	Doub a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
	Doub rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
	Int i, n = xold.size();
	check = false;
	for (i = 0; i<n; i++) sum += p[i] * p[i];
	sum = sqrt(sum);
	if (sum > stpmax)
		for (i = 0; i<n; i++)
			p[i] *= stpmax / sum;
	for (i = 0; i<n; i++)
		slope += g[i] * p[i];
	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	test = 0.0;
	for (i = 0; i<n; i++) {
		temp = abs(p[i]) / MAX(abs(xold[i]), 1.0);
		if (temp > test) test = temp;
	}
	alamin = TOLX / test;
	alam = 1.0;
	for (;;) {
		for (i = 0; i<n; i++) x[i] = xold[i] + alam * p[i];
		f = func(x);
		if (alam < alamin) {
			for (i = 0; i<n; i++) x[i] = xold[i];
			check = true;
			return;
		}
		else if (f <= fold + ALF * alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope / (2.0*(f - fold - slope));
			else {
				rhs1 = f - fold - alam * slope;
				rhs2 = f2 - fold - alam2 * slope;
				a = (rhs1 / (alam*alam) - rhs2 / (alam2*alam2)) / (alam - alam2);
				b = (-alam2 * rhs1 / (alam*alam) + alam * rhs2 / (alam2*alam2)) / (alam - alam2);
				if (a == 0.0) tmplam = -slope / (2.0*b);
				else {
					disc = b * b - 3.0*a*slope;
					if (disc < 0.0) tmplam = 0.5*alam;
					else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0*a);
					else tmplam = -slope / (b + sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam = 0.5*alam;
			}
		}
		alam2 = alam;
		f2 = f;
		alam = MAX(tmplam, 0.1*alam);
	}
}

template <class T>
struct NRfdjac {
	const Doub EPS;
	T &func;
	NRfdjac(T &funcc) : EPS(1.0e-8), func(funcc) {}
	MatDoub operator() (VecDoub_I &x, VecDoub_I &fvec) {
		Int n = x.size();
		MatDoub df(n, n);
		VecDoub xh = x;
		for (Int j = 0; j<n; j++) {
			Doub temp = xh[j];
			Doub h = EPS * abs(temp);
			if (h == 0.0) h = EPS;
			xh[j] = temp + h;
			h = xh[j] - temp;
			VecDoub f = func(xh);
			xh[j] = temp;
			for (Int i = 0; i<n; i++)
				df[i][j] = (f[i] - fvec[i]) / h;
		}
		return df;
	}
};

template <class T>
struct NRfmin {
	VecDoub fvec;
	T &func;
	Int n;
	NRfmin(T &funcc) : func(funcc) { n = 0; }
	Doub operator() (VecDoub_I &x) {
		n = x.size();
		Doub sum = 0;
		fvec = func(x);
		for (Int i = 0; i<n; i++) sum += SQR(fvec[i]);
		return 0.5*sum;
	}
};
double norm2(VecDoub x)
{
	double sqrtprod = 0;
	for (int i = 0; i< x.size(); i++)
	{
		sqrtprod += x[i] * x[i];
	}
	return sqrt(sqrtprod);
}

VecDoub diff(VecDoub x, VecDoub y)
{
	VecDoub res = VecDoub(x.size());
	for (int i = 0; i< x.size(); i++)
	{
		res[i] = x[i] - y[i];
	}
	return res;
}

template <class T>
void newt(VecDoub_IO &x, Bool &check, T &vecfunc) {
	const Int MAXITS = 2000;
	const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
	const Doub TOLX = numeric_limits<Doub>::epsilon();
	Int i, j, its, n = x.size();
	Doub den, f, fold, stpmax, sum, temp, test;
	VecDoub g(n), p(n), xold(n);
	MatDoub fjac(n, n);
	NRfmin<T> fmin(vecfunc);
	NRfdjac<T> fdjac(vecfunc);
	VecDoub &fvec = fmin.fvec;
	f = fmin(x);
	test = 0.0;
	for (i = 0; i<n; i++)
		if (abs(fvec[i]) > test) test = abs(fvec[i]);
	if (test < 0.01*TOLF) {
		check = false;
		return;
	}
	sum = 0.0;
	for (i = 0; i<n; i++) sum += SQR(x[i]);
	stpmax = STPMX * MAX(sqrt(sum), Doub(n));
	double dkOld = 0;
	for (its = 0; its<MAXITS; its++) {
		fjac = fdjac(x, fvec);
		for (i = 0; i<n; i++) {
			sum = 0.0;
			for (j = 0; j<n; j++) sum += fjac[j][i] * fvec[j];
			g[i] = sum;
		}
		cout << setw(15) << its + 1;
		for (j = 0; j<8; j++) cout << setw(15) << x[j];
		// cout << setw(15) << x[0];
		cout << setw(15) << norm2(diff(x, xold));
		cout << setw(15) << norm2(diff(x, xold)) / (dkOld*dkOld);
		cout << setw(15) << 'e';
		cout << endl;

		dkOld = norm2(diff(x, xold));


		for (i = 0; i<n; i++) xold[i] = x[i];
		fold = f;
		for (i = 0; i<n; i++) p[i] = -fvec[i];
		LUdcmp alu(fjac);
		alu.solve(p, p);
		lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);

		test = 0.0;
		for (i = 0; i<n; i++)
			if (abs(fvec[i]) > test) test = abs(fvec[i]);
		if (test < TOLF) {
			check = false;
			cout << "TOLF" << endl;
			return;
		}
		if (check) {
			test = 0.0;
			den = MAX(f, 0.5*n);
			for (i = 0; i<n; i++) {
				temp = abs(g[i])*MAX(abs(x[i]), 1.0) / den;
				if (temp > test) test = temp;
			}
			check = (test < TOLMIN);
			cout << "TOLMIN" << endl;
			return;
		}
		test = 0.0;
		for (i = 0; i<n; i++) {
			temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
			if (temp > test) test = temp;
		}
		if (test < TOLX)
			return;
	}
	throw("MAXITS exceeded in newt");
}



const double v = 120, k = 2.5, w = 4.0, a = 2e-7, d = 30.0;
double n;

VecDoub f(VecDoub_IO x)
{
	VecDoub res(8);

	res[0] = x[6] * (cosh(x[3] / x[6]) - 1) - x[2];  		// L0
	res[1] = 2 * x[6] * sinh(x[3] / x[6]) - x[1];			// L
	res[2] = 2 * x[3] + 2 * k*cos(x[4]) - d;			// p
	res[3] = x[2] + k * sin(x[4]) - n;			// x
	res[4] = sinh(x[3] / x[6]) - tan(x[5]);	// theta
	res[5] = (1 + v / (w*x[0]))*tan(x[5]) - tan(x[4]);	// phi
	res[6] = x[0] * (1 + a * x[7]) - x[1];			// a
	res[7] = w * x[0] / (2 * sin(x[5])) - x[7];			// H
	return res;
}

int main() {
	//	d = 30.0;
	vector<double> vec_n = { 5,2,1,0.5,0.2,0.1 };
	int j;
	bool check;
	VecDoub x(8);

	for (auto& i : vec_n)
	{
		cout << endl << "-----------------------------------------------------------------------------" << endl;
		cout << setw(15) << "k";
		for (int j = 0; j<8; j++) cout << setw(15) << "x_k";
		//cout << setw(15) << "x_k";
		cout << setw(15) << "dx_k";
		cout << setw(15) << "dx_k / dx_(k-1)²";
		cout << setw(15) << "e";
		cout << endl;

		n = i;

		x[0] = 30;			// L0
		x[1] = 30;			// L
		x[2] = 0.1;			// p
		x[3] = 15;			// x
		x[4] = 3.14 / 5;		// theta
		x[5] = 3.14 / 20;	// phi
		x[6] = 40;			// a
		x[7] = 5;			// H
							/*
							x[0] = 10;			// L0
							x[1] = 10;			// L
							x[2] = 10;			// p
							x[3] = 10;			// x
							x[4] = 1;		// theta
							x[5] = 1;		// phi
							x[6] = 10;			// a
							x[7] = 10;			// H

							*/


		newt(x, check, f);
		cout << "For n = " << i << " the result for L_0 is " << x[0] << endl;
	}

	cin.ignore();
	return 0;
}