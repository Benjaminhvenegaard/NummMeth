#include <iostream>


#include "../Source Code/code/nr3.h"
#include "../Source Code/code/utilities.h"
#include "../Source Code/code/banded.h"
#include "../Source Code/code/ludcmp.h"
using namespace std;

int kk = 1;
int n = pow(2, kk) + 1;
int nn = n * n;
Doub h = 1.0 / (n + 1);

MatDoub A(nn, nn);
MatDoub CompactA(nn, 2 * n + 1);
VecDoub b(nn);
VecDoub inx(nn);
VecDoub iny(nn);
VecDoub xx(nn);
VecDoub yy(nn);
VecDoub zz(nn);
Doub mid = (nn + 1) / 2 - 1;



int main() {
	std::cout << " n= " << n << endl;
	// Calculating indexes for matrix
	for (Int i = 0; i<nn; i++) {
		inx[i] = i % n + 1;
		iny[i] = floor(i / n) + 1;
	}
	util::print(inx, "inx");
	util::print(iny, "iny");
	// Calculating values for matrix
	for (Int i = 0; i<nn; i++) {
		xx[i] = inx[i] / (n + 1);
		yy[i] = iny[i] / (n + 1);
	}
	util::print(xx, "xx");
	util::print(yy, "yy");
	//Setting up A matrix and b vector
	for (Int i = 0; i<nn; i++) {
		A[i][i] = 4;
		int ix = inx[i];
		int iy = iny[i];

		if (ix > 1)
			A[i][i - 1] = -1.0;
		if (ix < n)
			A[i][i + 1] = -1.0;
		if (iy > 1)
			A[i][i - n] = -1.0;
		if (iy < n)
			A[i][i + n] = -1.0;

		b[i] = h * h*(1 + xx[i] + yy[i]);
	}

	util::print(A, "A");
	util::print(b, "b");
	//	Bandec ElipBand(A,n,n);
	//	ElipBand.solve(b,zz);
	//	util::print(zz,"Band zz");
	//	util::print(ElipBand.au,"au");
	//	util::print(ElipBand.al,"al");

	LUdcmp LU(A);
	LU.solve(b, zz);
	util::print(zz, "LU zz");

	cout << "zz mid" << endl;
	cout << zz[mid] << endl;

	cin.ignore();
	return 0;
}