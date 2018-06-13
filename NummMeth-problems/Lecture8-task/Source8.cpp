#include <iostream>


#include "../Source Code/code/nr3.h"
#include "../Source Code/code/quadrature.h"
#include "../Source Code/code/derule.h"


using namespace std;

double F00(double x)
{
	return cos(pow(x, 2))*exp(-x);
}

double F01(double x, double tau)
{
	return cos(pow(x, 2))*exp(-x);
}

double F02(double x, double tau)
{
	if (x < 0.01)
	{
		return (1 / sqrt(tau))*cos(pow(x, 2))*exp(-x);
	}
	return (1 / sqrt(x)) * cos(pow(x, 2))*exp(-x);
}

template<class T>
Doub trap(T &func, const Doub a, const Doub b, const Doub eps = 1.0e-10) {
	const Int JMAX = 20;
	Doub s, olds = 0.0, oldolds = 0.0, rn = 0.0, r = 0.0, oldr = 0.0, oldoldr = 0.0;
	Trapzd<T> t(func, a, b);
	for (Int j = 0; j<JMAX; j++) {
		s = t.next();
		cout << setw(15) << j + 1;
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


		if (j > 5)
			if (abs(s - olds) < eps*abs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;

		oldoldr = oldr;
		oldolds = olds;
		oldr = r;
		olds = s;

	}
	throw("Too many steps in routine qtrap");
}

template<class T>
double derule(T &func, const double a, const double b, const double pre)
{
	Doub count = 0;
	Doub N = 0;
	Doub Ah0 = 0;
	Doub Ah1 = 0;
	DErule<T> rule(func, a, b, 4.3); //3.7    4.3  Vælges ud fra sværheden af singulariteten ved 1/sqrt(x) er 4.3 passende
	while ((abs(Ah0 - Ah1) / Ah0 > pre || count == 0))//&& count < 16384 )
	{
		count++;
		N = pow(2, count);
		Ah1 = Ah0;
		Ah0 = rule.next();
		cout << count << "  " << Ah0 << "   " << abs(Ah0 - Ah1) << endl;
	}
	return Ah0;
}



int main() {
	cout << "cos(pow(x,2))*exp(-x)" << endl;
	cout << setw(15) << "n" << setw(15) << "Sn" << setw(15) << "k1" << setw(15) << "Rn" << setw(15) << "k2" << setw(15) << "R2n" << endl;
	cout << "result:" << trap(F00, 0, 1, 1e-10) << endl << endl;

	cout << "cos(pow(x,2))*exp(-x)" << endl;
	cout << "result:" << derule(F01, 0, 1, 1e-15) << endl << endl;

	cout << "(1/sqrt(tau))*cos(pow(x,2))*exp(-x) near x == 0" << endl;
	cout << "result:" << derule(F02, 0, 1, 1e-15) << endl;


	cin.ignore();
	return 0;
}