#include <iostream>

#include "../Source Code/code/nr3.h"
#include "../Source Code/code/tridag.h"
#include "../Source Code/code/utilities.h"
#include "psplot.h

struct xi {
	Doub a, h;
	xi(Doub aa, Doub hh) :a(aa), h(hh) {}
	Doub operator() (Doub i)
	{
		return a + i * h;
	}
};

struct yi {
	Doub a, b, N;
	yi(Doub aa, Doub bb, Doub NN) :a(aa), b(bb), N(NN) {}
	Doub operator() (Doub i)
	{
		return a + i / N * (b - a);
	}
};

double F(Doub x, Doub y, Doub yp)
{
	return -cos(y)*sin(yp);
	//return 2.0;
}
double Fy(Doub x, Doub y, Doub yp)
{
	return sin(y)*sin(yp);
	//return 0;
}
double Fyp(Doub x, Doub y, Doub yp)
{
	return -cos(y)*cos(yp);
	//return 0;
}




using namespace std;

int main() {
	// Initialize System
	Doub N = 100;
	Doub a = 0, b = 10, aa = 0, bb = 3;
	Doub h = (b - a) / N;
	xi x(a, h);
	yi yinit(aa, bb, N);
	VecDoub y(N + 1, 0.0);
	VecDoub dy(N + 1, 0.0);
	VecDoub Jp(N + 1, 0.0); //Note that Jp is one shorter than J
	VecDoub J(N + 1, 0.0);
	VecDoub Jm(N + 1, 0.0); //Note that Jm is one shorter than J
	VecDoub Fi(N + 1, 0.0);
	VecDoub Fm(N + 1, 0.0);
	for (int i = 0; i<N + 1; i++) { y[i] = yinit(i); } //Startguess for Y
	util::print(y, "y");

	//LOOP STARTS HERE
	for (int k = 0; k<4; k++) {
		//Update
		//N=10;
		//h=(b-a)/N;
		//x =xi(a,h);

		// Define Fi vector
		Fi[0] = y[0] - aa;
		Fi[N] = y[N] - bb;
		for (int i = 1; i<N; i++) { Fi[i] = y[i + 1] - 2 * y[i] + y[i - 1] - h * h*F(x(i), y[i], (y[i + 1] - y[i - 1]) / (h*2.0)); }
		for (int i = 0; i<N + 1; i++) { Fm[i] = -Fi[i]; }
		// Define JVectors
		J[0] = 1;
		J[N] = 1;
		for (int i = 1; i<N; i++)
		{
			Jm[i] = 1 + h / 2.0*Fyp(x(i), y[i], (y[i + 1] - y[i - 1]) / (h*2.0));        // Define the Ji,i-1 vector, called Jminus or Jm
			J[i] = -2 - h * h*   Fy(x(i), y[i], (y[i + 1] - y[i - 1]) / (h*2.0));        // Define the Ji,i vector, called J
			Jp[i] = 1 - h / 2.0*Fyp(x(i), y[i], (y[i + 1] - y[i - 1]) / (h*2.0));        // Define the Ji,i+1 vector, called Jplus or Jp
		}
		// if you want, put them in a Matrix.. alternatively use Tridag

		// util::print(Fi,"Fi");
		//util::print(Jp,"Jp");
		// util::print(J,"J");
		// util::print(Jm,"Jm");
		tridag(Jm, J, Jp, Fm, dy);
		for (int i = 0; i<N + 1; i++) { y[i] = y[i] + dy[i]; }
		util::print(y, "y");
	}
	// util::print(y,"y");

	PSpage pg("myplot.ps");
	PSplot plot1(pg, 100., 500., 500., 800.);
	//PSplot plot2(pg,100.,500.,100.,400.);
	plot1.setlimits(a, b, 0., 5.);
	plot1.frame();
	plot1.autoscales();
	plot1.xlabel("x");
	plot1.ylabel("y");
	//plot1.lineplot(xlag,ylag);

	for (int i = 0; i< y.size(); i++)
	{
		plot1.pointsymbol(x(i), y[i], 108, 8);
	}



	pg.close();
	pg.display("\"ghostscript\" ");
}