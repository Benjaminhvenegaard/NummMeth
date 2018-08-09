

#include<iostream>
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/ludcmp.h"
#include "../Source Code/code/utilities.h"


using namespace std;
using namespace util;

//Task for lecture 1


int main()
{
	MatDoub A(3, 3);

	A[0][0] = 1.0; A[0][1] = 2; A[0][2] = 3;
	A[1][0] = 2.0; A[1][1] = -4; A[1][2] = 6;
	A[2][0] = 3.0; A[2][1] = -9; A[2][2] = -3;

	VecDoub b(3);

	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	VecDoub x(3);

	print(A, "A");
	print(b, "b");
	print(x, "x");


	LUdcmp LU(A);

	auto L = LU.lu;

	L[2][1] = L[1][0] = L[2][0] = 0;
	L[0][0]   = L[1][1]  = L[2][2] = 1;

	print(L, "L");

	auto U = LU.lu;

	U[0][1] = U[0][2] = U[1][2] = 0;
	print(U, "U");



	// Solve LU object with b as result and put result in x
	LU.solve(b, x);

	print(A, "A");
	print(L*U, "L*U");
	print(L*U*x, "L*U*x");
	print(x, "x");


	cin.ignore();
	return 0;
}








/*
// NUM recipies includes:
#include "../Source Code/code/nr3.h"
#include "../Source Code/code/ludcmp.h"

// My own includes
#include "../Source Code/code/utilities.h"


using namespace std;
using namespace util;

// Exercise 1
// A x = b using LU decomposition

int main()
{
	// Make an A matrix which is a 3 by 3 matrix 
	MatDoub A(3, 3);
	// Insert values in the matrix
	A[0][0] = 1.0;	A[0][1] = 2.0;	A[0][2] = 3.0;
	A[1][0] = 2.0;	A[1][1] = -4.0;	A[1][2] = 6.0;
	A[2][0] = 3.0;	A[2][1] = -9.0;	A[2][2] = -3.0;

	// Make the b (result vector)
	VecDoub b(3);
	// Insert values in b-vector
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	// Make a variable-vector x
	VecDoub x(3);

	// Make the LU object with the matrix A naming it LU
	LUdcmp LU(A);

	// Solve the system

	auto L = LU.lu;
	L[0][1] = L[0][2] = L[1][2] = 0;
	L[0][0] = L[1][1] = L[2][2] = 1;
	print(L, "L");

	auto U = LU.lu;
	U[1][0] = U[2][0] = U[2][1] = 0;
	print(U, "U");

	// Solve LU object with b as result and put result in x
	LU.solve(b, x);

	print(A, "A");
	print(L*U, "L*U");
	print(L*U*x, "L*U*x");
	print(x, "x");





	cin.ignore();
	return 0;



}
*/
