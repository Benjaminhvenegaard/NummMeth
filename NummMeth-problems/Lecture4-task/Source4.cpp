
//Includes:
#include <iostream>
#include <fstream>
#include <math.h>


//H-files Numerical recipies:
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/svd.h"
#include "../Source Code/code/utilities.h"


using namespace std;
using namespace util;


int main()
{
	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("../Lecture2-task/FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	MatDoub A(xFilip.size(), 11);
	for (int i = 0; i <A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
			A[i][j] = pow(xFilip[i], (double)j);
	}


	SVD SVD_A = A;
	
	VecDoub x(A.ncols());
	VecDoub x2(A.ncols());
	
	
	MatDoub U = SVD_A.u;
	MatDoub V = SVD_A.v;
	VecDoub W = SVD_A.w;

	MatDoub Winv = MatDoub(W.size(), W.size());

	for (size_t i = 0; i < W.size(); i++)
	{
		for (size_t j = 0; j < W.size(); j++)
		{
			Winv[i][j] = 0;
		}
	}
	for (size_t i = 0; i < W.size(); i++)
	{
		Winv[i][i] = 1.0 / W[i];
	}

	

	x = V * Winv * (Transpose(U) * yFilip);

	SVD_A.solve(yFilip, x2, SVD_A.eps);
	cout << "Solution Filip using SVD:" << endl;
	print(x,"x with my own SVD");
	print(x2,"x with build in SVD");
	cout << endl;
	print(Winv, "Winv");



	cin.ignore();
	return 0;
}