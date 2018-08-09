
//Includes:
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

//H-files Numerical recipies:
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/ludcmp.h"
#include "../Source Code/code/utilities.h"
#include "../Source Code/code/cholesky.h"

//H-files My own:

//Namespaces
using namespace std;
using namespace util;

int main()
{

	VecDoub xPont(40); VecDoub yPont(40);
	ifstream Pont("../Lecture2-task/PontiusData.dat");
	for (int i = 0; i < 40; i++) 
	{
		Pont >> yPont[i];
		Pont >> xPont[i];
	}

	MatDoub A( 40, 3);

	for (int i = 0; i < A.nrows(); i++)
	{
		A[i][0] = 1.0;
		A[i][1] = xPont[i];
		A[i][2] = xPont[i] * xPont[i];
	}


	// AT* A * x = AT * y

	MatDoub At = Transpose(A);

	MatDoub AtA = At * A;
	VecDoub AtY = At * yPont;

	//Use LU to decompose

	LUdcmp LU(AtA);
	//Make the result vector
	VecDoub x(3);

	MatDoub C;

	LU.solve(AtY, x);

	print(x, "Pontius solution with LU");

	LU.inverse(C);

	//printDiag(C, "Variance");

	cout << endl;

	//Solve using Cholesky

	Cholesky Cho(AtA);

	Cho.solve(AtY, x);

	print(x, " Pontius solution using Cholesky:");

	Cho.inverse(C);
	printDiag(C, "Variance:");

	cout << endl << endl;

	cout << endl;
	cout << endl;





	//.....Filip-dataset.....//
	
	
	/*VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("../Lecture2-task/FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	MatDoub A2(xFilip.size(), 11);

	for (size_t i = 0; i < A2.nrows(); i++)
	{
		for (size_t j = 0; j < A2.ncols(); j++)
		{
			A2[i][j] = pow(xFilip[i], (double)j);
		}
	}



	// AT* A * x = AT * y

	MatDoub At2 = Transpose(A2);

	MatDoub AtA2 = At2 * A2;
	VecDoub AtY2 = At2 * yFilip;

	//Use LU to decompose

	VecDoub x2(A2.ncols());
	LUdcmp LU2(AtA2);
	//Make the result vector
	LU2.solve(AtY2, x2);

	print(x2, "Filip solution with LU");
	LU2.inverse(C);





	//printDiag(C, "Variance");

	cout << endl;

	//Solve using Cholesky

	Cholesky Cho2(AtA2);

	Cho2.solve(AtY2, x2);

	print(x2, " Pontius solution using Cholesky:");

	Cho2.inverse(C);
	printDiag(C, "Variance:");


	cout << endl << endl;
	*/
	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("../Lecture2-task/FilipData.dat");
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

	print(A2, "A2");
	MatDoub AT2 = util::Transpose(A2);

	// At * A * x = At * y
	MatDoub ATA2 = AT2 * A2;
	VecDoub ATy2 = AT2 * yFilip;

	VecDoub x2(A2.ncols());
	LUdcmp LUFilip(ATA2);
	LUFilip.solve(ATy2, x2);
	
	util::print(x2, "Filip using LU:");
	LUFilip.inverse(C);
		util::printDiag(C, "Variance:");
	cout << endl;
	cout << endl;
	cout << endl;






	cin.ignore();
	
}


















