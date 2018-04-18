#include <iostream>
#include <fstream>
#include "../Source Code/code/utilities.h"
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/quadrature.h"
#include "../Source Code/code/svd.h"
#include <math.h>

#define PI 3.14159265359
#define E 2.71828182845904523536
#define d 1.0
#define T1 1000
#define T2 500
#define epsilon1 0.80
#define epsilon2 0.60
#define rho 1.712*pow(10,-9)
#define w 1.0

using namespace std;
using namespace util;

Doub operator*(const VecDoub &a, const VecDoub &b)
{
	Doub res = 0;
	if (a.size() != b.size()) {
		cerr << "in prod: the number of rows in A is not equal to the size of vector b" << endl;
	}
	for (int i = 0; i < a.size(); i++)
	{
		res += a[i] * b[i];
	}
	return res;
}

void MatfindFpoints(Doub a, Doub b, Doub n, MatDoub &Matrix)
{
	Doub x, y;
	Doub value = 0;
	Doub h = (b-a) / (n);
	int nn = n + 1;
	vector<Doub> temp;
	MatDoub tempMatrix(nn, nn);
	for (int i = 0; i < nn; i++)
	{
		x = a + Doub(i)*h;
		for (int j = 0; j < nn; j++)
		{
			y = a + Doub(j)*h;
			value = 0.5 * (1.0 / (pow(pow(d, 2.0) + pow(x - y, 2.0), (3.0 / 2.0))));
			tempMatrix[i][j] = value;
		}
	}
	Matrix = tempMatrix;
}

void makeAMatrix(MatDoub &mat, MatDoub &fMat, Doub a, Doub b)
{
	int n = fMat.ncols();
	Doub h = (b - a) / (n-1);
	MatDoub tempMatrix(n * 2, n * 2);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (j+n == n || j+n == 2*n - 1)
				tempMatrix[i][j + n] = ((1.0 - epsilon1)*fMat[i][j])*h*0.5;
			else
				tempMatrix[i][j + n] = ((1.0 - epsilon1)*fMat[i][j]) * h;

			if (j == i)
				tempMatrix[i][j] = -1;
			else
				tempMatrix[i][j] = 0;
		}

	for (int i = n; i < 2 * n; i++)
		for (int j = 0; j < n; j++)
		{
			if (j == 0 || j  == n - 1)
				tempMatrix[i][j] = ((1.0 - epsilon2)*fMat[i - n][j]) * h*0.5;
			else
				tempMatrix[i][j] = ((1.0 - epsilon2)*fMat[i - n][j]) * h;
			
			
			if (j == i - n)
				tempMatrix[i][j+n] = -1;
			else
				tempMatrix[i][j+n] = 0;
		}

	mat = tempMatrix;
}

void makebVector(VecDoub &vec, MatDoub &mat)
{
	int n = mat.ncols();
	VecDoub tempVec(n);
	for (int i = 0; i < n; i++)
	{
		if (i < n / 2)
			tempVec[i] = -epsilon1*rho*pow(T1, 4);
		else
			tempVec[i] = -epsilon2*rho*pow(T2, 4);
	}
	vec = tempVec;
}

VecDoub trapezV(MatDoub &fMat, Doub a, Doub b, VecDoub &resVec)
{
	int n = fMat.ncols();
	Doub h = (b - a) / (n-1);
	VecDoub tempVec(n);
	Doub tempRes = 0;
	VecDoub res(n);
	VecDoub vec(n);
	for (int i = 0; i < n; i++)
		vec[i] = resVec[i + n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			tempVec[j] = fMat[i][j]; //ROW

		for (int j = 0; j < n; j++)
		{
			if (j == 0)
				tempRes += 0.5*(tempVec[j]*vec[i])*h;
			else if (j == (n - 1))
				tempRes += 0.5*(tempVec[j]*vec[i])*h;
			else
				tempRes += (tempVec[j]*vec[i])*h;
		}
		res[i] = tempRes;
		tempRes = 0;
	}
	return res;
}

VecDoub trapezU(MatDoub &fMat, Doub a, Doub b, VecDoub &resVec)
{
	int n = fMat.ncols();
	Doub h = (b - a) / (n-1);
	VecDoub tempVec(n);
	Doub tempRes = 0;
	VecDoub res(n);
	VecDoub vec(n);
	for (int i = 0; i < n; i++)
		vec[i] = resVec[i];
	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
			tempVec[j] = fMat[j][i]; //COLLUM

		for (int j = 0; j < n; j++)
		{
			if (j == 0)
				tempRes += 0.5*(tempVec[j] * vec[i])*h;
			else if (j == (n - 1))
				tempRes += 0.5*(tempVec[j] * vec[i])*h;

			else
				tempRes += (tempVec[j] * vec[i])*h;
		}
		res[i] = tempRes;
		tempRes = 0;
	}
	return res;
}

VecDoub trapezQ(Doub a, Doub b, VecDoub &I1Vec, VecDoub &I2Vec, VecDoub &resVec)
{
	int n = I1Vec.size();
	Doub h = (b - a) / (n-1);
	VecDoub tempVec(n);
	Doub tempRes = 0;
	Doub Q1 = 0;
	Doub Q2 = 0;
	VecDoub QVec(2);
	VecDoub uVec(n);
	VecDoub vVec(n);

	for (int i = 0; i < n; i++)
		uVec[i] = resVec[i];

	for (int i = 0; i < n; i++)
		vVec[i] = resVec[i+n];

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			Q1 += 0.5*(uVec[i] - I1Vec[i])*h;
			Q2 += 0.5*(vVec[i] - I2Vec[i])*h;
		}
		else if (i == (n - 1))
		{
			Q1 += 0.5*(uVec[i] - I1Vec[i])*h;
			Q2 += 0.5*(vVec[i] - I2Vec[i])*h;
		}

		else
		{
			Q1 += (uVec[i] - I1Vec[i])*h;
			Q2 += (vVec[i] - I2Vec[i])*h;
		}
	}
	QVec[0] = Q1;
	QVec[1] = Q2;
	return QVec;
}


int main() {
	VecDoub I1, I2, QVec;
	vector<vector<Doub>> fVec;   //Indeholder alle funktionsværdier af F
	MatDoub fMatrix;			//
	MatDoub AMatrix;			//
	VecDoub bVec;				
	Doub a = - 0.5;				//Endepunkterne
	Doub b = 0.5;				//
	Doub alpha = 2, k = 2;		//
	Doub s1=0, s2=0, lastS1=0, lastS2=0, prevLastS1=0, prevLastS2=0;

	for (int n = 4; n <= 128; n = n * 2)
	{
		prevLastS1 = lastS1;
		prevLastS2 = lastS2;
		lastS1 = s1;
		lastS2 = s2;

		MatfindFpoints(-0.5, 0.5, n, fMatrix);
		makeAMatrix(AMatrix, fMatrix, -0.5, 0.5);
		
		makebVector(bVec, AMatrix);
		
		VecDoub resVec(bVec.size());
		
		SVD svd(AMatrix);
		svd.solve(bVec, resVec, svd.eps);
		//util::print(resVec, "Result");
		//print(AMatrix);
		//print(bVec);
		//print(fMatrix);
		I1 = trapezV(fMatrix, a, b, resVec);
		I2 = trapezU(fMatrix, a, b, resVec);
		QVec = trapezQ(a, b, I1, I2, resVec);
		cout << endl;
		cout << setw(5) << "u(-0.5): " << resVec[0] << setw(15) << "u(-0.25): "  << resVec[(resVec.size() - 1)*0.125] << setw(15) << "u(0): " << resVec[(resVec.size() - 1) / 4] << setw(15) << "u(0.25): "  << resVec[(resVec.size() - 1)*0.375] << setw(15) << "u(0.5): " << resVec[(resVec.size() - 1)/2] << endl;
		cout << setw(5) << "v(-0.5): " << resVec[resVec.size()/2] << setw(15) << "v(-0.25): " << resVec[(resVec.size())*0.625] << setw(15) << "v(0): " << resVec[(resVec.size()) *0.75] << setw(15) << "v(0.25): " << resVec[(resVec.size())*0.875] << setw(15) << "v(0.5): " << resVec[(resVec.size()-1)] << endl;
		cout << endl;
		cout << setw(5) << "n" << setw(15) << "Q1" << setw(15) << "Q2" << setw(15) << "error1" << setw(15) << "error2" << setw(15) << "Order1" << setw(15) << "Order2" << endl;
		//cout << "Q for " << n << ": " << endl;
		//print(QVec, "QVec");
		s1 = QVec[0];
		s2 = QVec[1];
		cout << setw(5) << n << setw(15) << QVec[0] << setw(15) << QVec[1] << setw(15) << (s1 - lastS1) / (pow(alpha, k) - 1) << setw(15) << (s2 - lastS2) / (pow(alpha, k) - 1) << setw(15) << (prevLastS1 - lastS1) / (lastS1 - s1) << setw(15) << (prevLastS2 - lastS2) / (lastS2 - s2) << endl;
		cout << endl;
		cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		//cout << setw(5) << N << setw(15) << s << setw(15) << (s - lastS) / (pow(alpha, k) - 1) << setw(15) << (lastLastS - lastS) / (lastS - s) << endl;
	}
	/*
	cout << "I1: " << endl;
	print(I1,"I1");
	cout << "I2: " << endl;
	print(I2, "I2");
	*/
	std::cout << "done" << endl;
	cin.ignore();
	return 0;
}