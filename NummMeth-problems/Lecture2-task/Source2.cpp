
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include "../Source Code/code/nr3.h"
#include "../Source Code/code/ludcmp.h"
#include "../Source Code/code/cholesky.h"
#include "../Source Code/code/svd.h"
#include "../Source Code/code/utilities.h"

using namespace std;

int main() {

	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("src/FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	// Get information from the data file
	VecDoub xPont(40); VecDoub yPont(40);
	ifstream Pont("src/PontiusData.dat");
	for (int i = 0; i < 40; i++) {
		Pont >> yPont[i];
		Pont >> xPont[i];
	}

	// your code
	// Make an A matrix with the informartion from the data file
	MatDoub A(40, 3);
	for (int i = 0; i <A.nrows(); i++)
	{
		A[i][0] = 1.0; A[i][1] = xPont[i]; A[i][2] = xPont[i] * xPont[i];
	}

	// Make the transpose of A
	MatDoub AT = util::Transpose(A);

	// At * A * x = At * y
	MatDoub ATA = AT * A;
	VecDoub ATy = AT * yPont;

	// Use LU decompose
	LUdcmp LU(ATA);
	// make a result vector
	VecDoub x(3);

	MatDoub C;
	LU.solve(ATy, x);

	util::print(x, "Pontius solution using LU:");
	LU.inverse(C);
	//	util::printDiag(C,"Variance:");
	cout << endl;

	// Solve using cholesky
	Cholesky Cho(ATA);
	Cho.solve(ATy, x);
	util::print(x, "Pontius solution using Cholesky:");

	Cho.inverse(C);
	//	util::printDiag(C,"Variance:");
	cout << endl << endl;
	/*
	* Cholesky
	*/

	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("src/FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	MatDoub A2(xFilip.size(), 11);
	for (int i = 0; i <A2.nrows(); i++)
	{
		for (int j = 0; j < A2.ncols(); j++)
			A2[i][j] = pow(xFilip[i], (double)j);
	}

	MatDoub AT2 = util::Transpose(A2);

	// At * A * x = At * y
	MatDoub ATA2 = AT2 * A2;
	VecDoub ATy2 = AT2 * yFilip;

	VecDoub x2(A2.ncols());
	LUdcmp LUFilip(ATA2);
	LUFilip.solve(ATy2, x2);
	util::print(x2, "Filip using LU:");
	LUFilip.inverse(C);
	//	util::printDiag(C, "Variance:");
	cout << endl;

	cout << "Filip solution using Cholesky:" << endl;
	try {
		Cholesky CholeskyFilip(ATA2);
	}
	catch (NRerror e) { cout << e.message << endl; }





	return 0;
}