

//Includes:
#include <iostream>
#include <math.h>


//H-files Numerical recipies:
#include "../Source Code/code/nr3.h"

using namespace std;

void line()
{
	cout << "-------------------------------------------------------------------------------------------" << endl;
}


struct f
{
	Doub operator()(Doub x) { return x - cos(x); }
	Doub dc(Doub x) { return 1 + sin(x); }

};


Doub pow10(int someInt)
{
	return pow(10,someInt);
}


// Make bisection method for root finding on non-linear equations
template <class T>
Doub Bisec(T &func, Doub x0, Doub x1, Doub acc)
{

}



int main()
{









	cin.ignore();
	return 0;
}


/*
#include <iostream>
#include "../Source Code/code/nr3.h"
#include <math.h>
//#include <roots.h>

using namespace std;

void line() { cout << endl << "-----------------------------------------------------------------------------------------------" << endl; }

//double f(double x) {return x-cos(x); }

struct f
{
Doub operator()(Doub x) { return x - cos(x); }
Doub dc(Doub x) { return 1 + sin(x); }
};

Doub pow10(int someint)
{
return (pow(10, someint));
}

// Make bisection method for root finding on non-linear equations
template <class T>
Doub rtbis(T &func, const Doub x1, const Doub x2, const Doub xacc)
{
const Int JMAX = 50;
Doub dx, xmid, rtb, dxOld = 99999;
Doub f = func(x1);
Doub fmid = func(x2);
if (f*fmid >= 0.0) throw("Root must be bracketed for bisection in rtbis");
rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
for (Int j = 0; j<JMAX; j++)
{
cout << setw(15) << j + 1;
cout << setw(15) << rtb;
cout << setw(15) << rtb + dx;
cout << setw(15) << dx;
cout << setw(15) << dx / dxOld;
cout << setw(15) << 'e' << endl;
dxOld = dx;


fmid = func(xmid = rtb + (dx *= 0.5));
if (fmid <= 0.0) rtb = xmid;
if (abs(dx) < xacc || fmid == 0.0) return rtb;
}
throw("Too many bisections in rtbis");
}


// Make secant method
template <class T>
Doub rtsec(T &func, const Doub x1, const Doub x2, const Doub xacc) {
const Int MAXIT = 30;
Doub xl, rts, xOld = 9999;
Doub fl = func(x1);
Doub f = func(x2);
if (abs(fl) < abs(f)) {
rts = x1;
xl = x2;
SWAP(fl, f);
}
else {
xl = x1;
rts = x2;
}
for (Int j = 0; j<MAXIT; j++) {
Doub dx = (xl - rts)*f / (f - fl);


cout << setw(15) << j + 1;
cout << setw(15) << rts;
cout << setw(15) << dx;
cout << setw(15) << abs(dx) / pow(abs(xOld - rts), 1.62);
cout << setw(15) << 'e' << endl;
xOld = rts;

xl = rts;
fl = f;
rts += dx;
f = func(rts);
if (abs(dx) < xacc || f == 0.0) return rts;
}
throw("Maximum number of iterations exceeded in rtsec");
}

template <class T>
Doub rtflsp(T &func, const Doub x1, const Doub x2, const Doub xacc) {
const Int MAXIT = 30;
Doub xl, xh, del, delOldl = 0, delOldh = 0;
Doub fl = func(x1);
Doub fh = func(x2);
if (fl*fh > 0.0) throw("Root must be bracketed in rtflsp");
if (fl < 0.0) {
xl = x1;
xh = x2;
}
else {
xl = x2;
xh = x1;
SWAP(fl, fh);
}
Doub dx = xh - xl;
for (Int j = 0; j<MAXIT; j++) {
Doub rtf = xl + dx * fl / (fl - fh);
Doub f = func(rtf);

if (f < 0.0) {
del = xl - rtf;
xl = rtf;
fl = f;

cout << setw(15) << j + 1;
cout << setw(15) << xl;
cout << setw(15) << xh;
cout << setw(15) << del;
cout << setw(15) << del / delOldl;
cout << setw(15) << 'e' << endl;
delOldl = del;
delOldh = 0;
}
else {
del = xh - rtf;
xh = rtf;
fh = f;

cout << setw(15) << j + 1;
cout << setw(15) << xl;
cout << setw(15) << xh;
cout << setw(15) << del;
cout << setw(15) << del / delOldh;
cout << setw(15) << 'e' << endl;
delOldl = 0;
delOldh = del;
}


dx = xh - xl;
if (abs(del) < xacc || f == 0.0) return rtf;
}
throw("Maximum number of iterations exceeded in rtflsp");
}

template <class T>
Doub zriddr(T &func, const Doub x1, const Doub x2, const Doub xacc) {
const int MAXIT = 60;
Doub fl = func(x1);
Doub fh = func(x2);
if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
Doub dx_old = 9.99e99;
Doub xl = x1;
Doub xh = x2;
Doub ans = -9.99e99;
for (Int j = 0; j<MAXIT; j++)
{
Doub xm = 0.5*(xl + xh);
Doub fm = func(xm);
Doub s = sqrt(fm*fm - fl * fh);
if (s == 0.0) return ans;
Doub xnew = xm + (xm - xl)*((fl >= fh ? 1.0 : -1.0)*fm / s);
if (abs(xnew - ans) <= xacc) return ans;

Doub dx = xnew - ans;
cout << setw(15) << j + 1;
cout << setw(15) << xnew;
cout << setw(15) << dx;
cout << setw(15) << abs(dx) / pow(abs(dx_old), 2);
cout << setw(15) << 'e' << endl;
dx_old = dx;

ans = xnew;
Doub fnew = func(ans);
if (fnew == 0.0) return ans;
if (SIGN(fm, fnew) != fm) {
xl = xm;
fl = fm;
xh = ans;
fh = fnew;
}
else if (SIGN(fl, fnew) != fl) {
xh = ans;
fh = fnew;
}
else if (SIGN(fh, fnew) != fh) {
xl = ans;
fl = fnew;
}
else throw("never get here.");


if (abs(xh - xl) <= xacc) return ans;
}
throw("zriddr exceed maximum iterations");
}
else {
if (fl == 0.0) return x1;
if (fh == 0.0) return x2;
throw("root must be bracketed in zriddr.");
}
}

template <class T>
Doub rtnewt(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
const Int JMAX = 20;
Doub rtn = 0.5*(x1 + x2), dxOld = 9999;
for (Int j = 0; j<JMAX; j++) {
Doub f = funcd(rtn);
Doub df = funcd.dc(rtn);
Doub dx = f / df;

cout << setw(15) << j + 1;
cout << setw(15) << rtn;
cout << setw(15) << dx;
cout << setw(15) << abs(dx) / pow(abs(dxOld), 2);
cout << setw(15) << 'e' << endl;
dxOld = dx;

rtn -= dx;
if ((x1 - rtn)*(rtn - x2) < 0.0)
throw("Jumped out of brackets in rtnewt");
if (abs(dx) < xacc) return rtn;
}
throw("Maximum number of iterations exceeded in rtnewt");
}


int main() {
//f function = f();
f function;
cout << "Bisection" << endl;
cout << setw(15) << "k" << setw(15) << "xmin" << setw(15) << "xmax" << setw(15) << "dx" << setw(15) << "C" << setw(15) << "e" << endl;
cout << "Result: " << rtbis(function, 0, 3.14 / 2, pow10(-8)) << endl;

line();
cout << "Secant" << endl;
cout << setw(15) << "k" << setw(15) << "x" << setw(15) << "dx" << setw(15) << "C" << setw(15) << "e" << endl;
cout << "Result: " << rtsec(function, 0, 3.14 / 2, pow10(-16)) << endl;

line();
cout << "False Position" << endl;
cout << setw(15) << "k" << setw(15) << "xmin" << setw(15) << "xmax" << setw(15) << "dx" << setw(15) << "C" << setw(15) << "e" << endl;
cout << "Result: " << rtflsp(function, 0, 3.14 / 2, pow10(-16)) << endl;

line();
cout << "Ridder" << endl;
cout << setw(15) << "k" << setw(15) << "x" << setw(15) << "dx" << setw(15) << "C" << setw(15) << "e" << endl;
cout << "Result: " << zriddr(function, 0, 3.14 / 2, pow10(-16)) << endl;

line();
cout << "Newton" << endl;
cout << setw(15) << "k" << setw(15) << "x" << setw(15) << "dx" << setw(15) << "C" << setw(15) << "e" << endl;
cout << "Result: " << rtnewt(function, -3.14, 3.14, pow10(-16)) << endl;



cin.ignore();
return 0;
}
*/