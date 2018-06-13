#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "../Source Code/code/nr3.h"


double Euler(double x0, double y0, int a, int b, int hinv)
{

	double x = x0;
	double y = y0;
	double t = a;
	double stepsize = (double(b) - double(a)) / double(hinv);
	//std::cout << "Euler method, with stepsize " << stepsize << std::endl;
	//std::cout << "x(0) = " << x << std::endl;
	//std::cout << "y(0) = " << y << std::endl;
	//std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
	//std::cout << std::endl;


	//for(double i = a; i< b-stepsize; i = i+stepsize)
	for (double i = a; i< b; i = i + stepsize)
	{
		double xplus1 = x + stepsize * (x*y);
		double yplus1 = y + stepsize * (-pow(x, 2));
		x = xplus1;
		y = yplus1;
		t = t + stepsize;
		//std::cout << "x(" << t << ") = " << x << std::endl;
		//std::cout << "y(" << t << ") = " << y << std::endl;
		//std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
		//std::cout << std::endl;

	}
	//std::cout << "x(" << t << ") = " << x << std::endl;
	return(x);
}

void Leap_Frog(double x0, double y0, int a, int b, int hinv)
{

	double x = x0;
	double y = y0;

	double t = a;
	double stepsize = (double(b) - double(a)) / double(hinv);

	double xminus1 = x - stepsize * (x*y);		// Using Euler
	double yminus1 = y - stepsize * (-pow(x, 2));	// Using Euler

	std::cout << "Leap_Frog method, with stepsize " << stepsize << std::endl;
	std::cout << "x(0) = " << x << std::endl;
	std::cout << "y(0) = " << y << std::endl;
	std::cout << "x^2 + y^2 = " << pow(x, 2) + pow(y, 2) << std::endl;
	std::cout << std::endl;

	for (double i = a; i< b; i = i + stepsize)
	{
		double xplus1 = xminus1 + 2 * stepsize*(x*y);
		double yplus1 = yminus1 + 2 * stepsize*(-pow(x, 2));
		xminus1 = x;
		yminus1 = y;
		x = xplus1;
		y = yplus1;
		t = t + stepsize;
		std::cout << "x(" << t << ") = " << x << std::endl;
		std::cout << "y(" << t << ") = " << y << std::endl;
		std::cout << "x^2 + y^2 = " << pow(x, 2) + pow(y, 2) << std::endl;
		std::cout << std::endl;

	}


}



double Midpoint(double x0, double y0, int a, int b, int hinv)
{

	double x = x0;
	double y = y0;

	double t = a;
	double stepsize = (double(b) - double(a)) / double(hinv);
	//double stepsize = 1.0/double(hinv);

	//std::cout << "Midpoint method, with stepsize " << stepsize << std::endl;
	//std::cout << "x(0) = " << x << std::endl;
	//std::cout << "y(0) = " << y << std::endl;
	//std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
	//std::cout << std::endl;


	for (double i = a; i< b; i = i + stepsize)
	{
		double xhalf = x + 0.5*stepsize*(x*y);
		double yhalf = y + 0.5*stepsize*(-pow(x, 2));

		double xplus1 = x + stepsize * (xhalf*yhalf);
		double yplus1 = y + stepsize * (-pow(xhalf, 2));
		x = xplus1;
		y = yplus1;
		t = t + stepsize;
		//std::cout << "x(" << t << ") = " << x << std::endl;
		//std::cout << "y(" << t << ") = " << y << std::endl;
		//std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
		//std::cout << std::endl;

	}
	return (x);

}


void Trapezoidal(double x0, double y0, int a, int b, int hinv)
{

	double x = x0;
	double y = y0;

	double t = a;
	double stepsize = (double(b) - double(a)) / double(hinv);

	std::cout << "Trapezoidal method, with stepsize " << stepsize << std::endl;
	std::cout << "x(0) = " << x << std::endl;
	std::cout << "y(0) = " << y << std::endl;
	std::cout << "x^2 + y^2 = " << pow(x, 2) + pow(y, 2) << std::endl;
	std::cout << std::endl;


	for (double i = a; i< b; i = i + stepsize)
	{
		double xpluseuler = x + stepsize * (x*y);
		double ypluseuler = y + stepsize * (-pow(x, 2));

		double xplus1 = x + 0.5*stepsize*((xpluseuler*ypluseuler) + (x*y));
		double yplus1 = y + 0.5*stepsize*((-pow(xpluseuler, 2)) + (-pow(x, 2)));
		x = xplus1;
		y = yplus1;
		t = t + stepsize;
		std::cout << "x(" << t << ") = " << x << std::endl;
		std::cout << "y(" << t << ") = " << y << std::endl;
		std::cout << "x^2 + y^2 = " << pow(x, 2) + pow(y, 2) << std::endl;
		std::cout << std::endl;

	}


}




int main() {
	Doub s, olds = 0.0, oldolds = 0.0, rn = 0.0, r = 0.0, oldr = 0.0, oldoldr = 0.0;
	cout << setw(15) << "hinv" << setw(15) << "Sn" << setw(15) << "k1" << setw(15) << "Rn" << setw(15) << "k2" << setw(15) << "R2n" << endl;
	for (int j = 0; j<10; j++) {


		//s = Euler(1,1,0,5,5*pow(2,j));
		//Leap_Frog(1,1,0,10,40);

		//s = Midpoint(1, 1, 0, 10, 5 * pow(2, j));

		Trapezoidal(1,1,0,10,40);

		cout << setw(15) << 5 * pow(2, j);
		/*
		cout << setw(15) << s;
		if (oldolds != 0.0) { cout << setw(15) << log2(abs((oldolds - olds) / (olds - s))); }
		else { cout << setw(15) << " "; }
		if (olds != 0.0) { r = s + (s - olds) / 3; cout << setw(15) << r; } // We see that k1 tends towards 2,and with alpha=2 we get alpha^k-1 = 3
		else { cout << setw(15) << " "; }
		if (oldoldr != 0.0) { cout << setw(15) << log2(abs((oldoldr - oldr) / (oldr - r))); }
		else { cout << setw(15) << " "; }
		if (oldr != 0.0) { rn = r + (r - oldr) / 15; cout << setw(15) << rn; } // We see that k2 tends towards 4,and with alpha=2 we get alpha^k-1 = 15
		else { cout << setw(15) << " "; }
		cout << endl;


		oldoldr = oldr;
		oldolds = olds;
		oldr = r;
		olds = s;
		*/

	}

	cin.ignore();
	return 0;
}