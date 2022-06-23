
//Lektion 1 opgave:
/*
#include <iostream>
//Lection 1 Includes:
#include "LUDecomp.h"

int main()
{
// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.
	//kan evt. lave en int variable size og sætte den til størrelsen af matrixen og vectoren.
	MatDoub A(3,3);
	A[0][0] = 1.0;	A[0][1] = 2.0;	A[0][2] = 3.0;
	A[1][0] = 2.0;	A[1][1] = -4.0;	A[1][2] = 6.0;
	A[2][0] = 3.0;	A[2][1] = -9.0;	A[2][2] = -3.0;

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	LUDecomp::LUDecomp(A, b, 3);
}*/



//Lektion 2 opgave:
/*
#include <iostream>
#include <fstream>
#include "LUDecomp.h"
#include "CholeskyDecomp.h"
#include "SVDDecomp.h"

using namespace std;

int main() {
	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("C:\\Users\\nvigg\\Documents\\GitHub\\Numeriskemetoder\\Data\\FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}
	std::cout << "Filip loaded" << std::endl;
	
	VecDoub xPont(40); VecDoub yPont(40);
	ifstream Pont("C:\\Users\\nvigg\\Documents\\GitHub\\Numeriskemetoder\\Data\\PontiusData.dat");
	for (int i = 0; i < 40; i++) {
		Pont >> yPont[i];
		Pont >> xPont[i];
	}
	std::cout << "Pontius loaded" << std::endl;
	
	// your code
	//Pontius
	MatDoub aPont(40, 3);
	double sigma = 0.205177424076185E-03;
	for (int i = 0; i < xPont.size(); ++i) {
		aPont[i][0] = 1. / sigma;
		aPont[i][1] = xPont[i] / sigma;
		aPont[i][2] = pow(xPont[i], 2) / sigma;
	}
	std::cout << "Matrix A made" << std::endl;
	VecDoub bPont(40);
	for (int i = 0; i < yPont.size(); ++i) {
		bPont[i] = yPont[i] / sigma;
	}
	std::cout << "Vector b made" << std::endl;

	//LU Decomposition
	MatDoub CPont;
	VecDoub cPont;

	CPont = util::Transpose(aPont) * aPont;
	cPont = util::Transpose(aPont) * bPont;

	LUDecomp::LUDecompquick(CPont, cPont, cPont.size());


	//Cholesky
	CholeskyDecomp::CholeskyDecomp(CPont, cPont, cPont.size());


	//Filip
	MatDoub aFilip(82, 11);
	double sigmaF = 0.334801051324544E-02;
	for (int i = 0; i < xFilip.size(); ++i) {
		aFilip[i][0] = 1 / sigmaF;
		aFilip[i][1] = xFilip[i] / sigmaF;
		aFilip[i][2] = pow(xFilip[i], 2) / sigmaF;
		aFilip[i][3] = pow(xFilip[i], 3) / sigmaF;
		aFilip[i][4] = pow(xFilip[i], 4) / sigmaF;
		aFilip[i][5] = pow(xFilip[i], 5) / sigmaF;
		aFilip[i][6] = pow(xFilip[i], 6) / sigmaF;
		aFilip[i][7] = pow(xFilip[i], 7) / sigmaF;
		aFilip[i][8] = pow(xFilip[i], 8) / sigmaF;
		aFilip[i][9] = pow(xFilip[i], 9) / sigmaF;
		aFilip[i][9] = pow(xFilip[i], 10) / sigmaF;
	}
	std::cout << "Matrix A made" << std::endl;

	VecDoub bFilip(82);
	for (int i = 0; i < yFilip.size(); ++i) {
		bFilip[i] = yFilip[i] / sigmaF;
	}
	std::cout << "Vektor B made" << std::endl;

	//LU Decomposition
	MatDoub CFilip;
	VecDoub cFilip;

	CFilip = util::Transpose(aFilip) * aFilip;
	cFilip = util::Transpose(aFilip) * bFilip;

	LUDecomp::LUDecompquick(CFilip, cFilip, cFilip.size());


	//Cholesky
	CholeskyDecomp::CholeskyDecomp(CFilip, cFilip, cFilip.size());
	
	//Dette deviere fra alle de andre svar af en eller anden grund
	
	//Lektion 3 SVD-decomposition af filip og pontius
	cout << "Pontius SVD:" << endl;
	SVDDecomp::SVDDecomp(aPont, bPont, aPont.ncols());

	cout << "Filip SVD:" << endl;
	SVDDecomp::SVDDecomp(aFilip, bFilip, aFilip.ncols());
	//Problemer med at få W til at indeholde noget!!!!
	SVDDecomp::SVDSolveAlg(aFilip, bFilip, aFilip.ncols());

	return 0;
}
*/



//Lektion 5 + Lektion 6
/*
#include "ConvergenceMethods.h"

using namespace std;

double f(double x) {
	return x - cos(x);
}

double df(double x) {
	return 1 - sin(x);
}


int main() {

	//Newtons method
	cout << "Newtons Method:" << endl;
	ConvergenceMethods::Newton(f, df, 10);

	//Secant method
	cout << "Secant Method:" << endl;
	ConvergenceMethods::Secant(f, df, 8);

	//Bitsection method
	cout << "Bitsection Method:" << endl;
	ConvergenceMethods::Bitsection(f, df, 10);

	//False position method
	cout << "False Position (Regula Falsi) Method:" << endl;
	ConvergenceMethods::RegulaFalsi(f, df, 10);

	//Ridders method
	cout << "Ridders Method:" << endl;
	ConvergenceMethods::Ridders(f, df, 5);

	return 0;
}
*/


//Lektion 7
/*
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "..\..\NR_C301\code\nr3.h"
#include "..\..\NR_C301\code\ludcmp.h"
#include "..\..\NR_C301\code\qrdcmp.h"
#include "..\..\NR_C301\code\roots_multidim.h"
#include "..\..\NR_C301\code\utilities.h"

double n = 5;

void printOutput(MatDoub& x, double n) {
	//x[0] = L0, x[1] = L, x[2] = p, x[3] = x, x[4] = thetha, x[5] = phi, x[6] = a, x[7] = H
	cout << left;
	cout << fixed;
	cout << setprecision(7);

	cout << setw(12) << "n" << " | " << setw(12) << "L0" << " | " << setw(12) << "L" << " | " <<
		setw(12) << "p" << " | " << setw(13) << "x" << " | " <<
		setw(13) << "thetha" << " | " << setw(14) << "phi" << " | " <<
		setw(14) << "a" << " | " << setw(12) << "H" << " | " << endl << endl;

		util::print(x);
}

VecDoub equations(VecDoub x) {
	//x[0] = L0, x[1] = L, x[2] = p, x[3] = x, x[4] = thetha, x[5] = phi, x[6] = a, x[7] = H
	double v = 120;
	double k = 2.5;
	double w = 4.0;
	double alpha = 2*10E-7;
	double d = 30;
	//double n = 5;
	

	VecDoub y(8);
	y[0] = x[6] * (cosh(x[3] / x[6]) - 1) - x[2];
	y[1] = 2 * x[6] * sinh(x[3] / x[6]) - x[1];
	y[2] = 2 * x[3] + 2 * k * cos(x[4]) - d;
	y[3] = x[2] + k * sin(x[4]) - n;
	y[4] = sinh(x[3] / x[6]) - tan(x[5]);
	y[5] = (1 + (v / (w * x[0]))) * tan(x[5]) - tan(x[4]);
	y[6] = x[0] * (1 + alpha * x[7]) - x[1];
	y[7] = ((w * x[0]) / (2 * sin(x[5]))) - x[7];

	return y;
}


int main() {
	VecDoub x(8);
	x[0] = 30;			//L0
	x[1] = 30;			//L
	x[2] = 0.1;			//P
	x[3] = 15;			//x
	x[4] = 3.14/5;		//Theta
	x[5] = 3.14/20;		//phi
	x[6] = 40;			//a
	x[7] = 5;			//H

	MatDoub j(5,9);
	bool convergence = false;
	for (int i = 0; i < 5; ++i) {
		n = 5 - i;
		newt(x, convergence, equations);
		j[i][0] = n;
		for (int k = 1; k < 9; ++k) {
			j[i][k] = x[k-1];
		}
	}
	printOutput(j, n);

	return 0;
}
*/

/*
//Lektion 8
#include "NewtonsCotes.h"

double f1(double x) {
	double ans = cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f2(double x) {
	double ans = sqrt(x) * cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f3(double x) {
	double ans = (1 / sqrt(x)) * cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f4(double x) {
	double ans = 1000 * exp(-1 / x) * exp(-1 / (1 - x));
	return ans;
}

int main() {
	VecDoub a(10);
	VecDoub b(10);
	VecDoub c(10);
	VecDoub d(10);
	a.assign(10, 0);
	b.assign(10, 0);
	c.assign(10, 0);
	d.assign(10, 0);
	VecDoub N(10);
	for (int i = 0; i < 10; ++i) {
		N[i] = pow(2, i) + 1;
		a[i] = NewtonsCotes::rectangleMethod(f1, N[i], 0, 1);
		b[i] = NewtonsCotes::rectangleMethod(f2, N[i], 0, 1);
		c[i] = NewtonsCotes::rectangleMethod(f3, N[i], 0, 1);
		d[i] = NewtonsCotes::rectangleMethod(f4, N[i], 0, 1);
	}
	std::cout << "Rectangle Method: " << std::endl;
	std::cout << "f1: " << std::endl;
	NewtonsCotes::printoutTable(a, N);
	std::cout << "f2: " << std::endl;
	NewtonsCotes::printoutTable(b, N);
	std::cout << "f3: " << std::endl;
	NewtonsCotes::printoutTable(c, N);
	std::cout << "f4: " << std::endl;
	NewtonsCotes::printoutTable(d, N);
	
	std::cout << "Trapezoidal Method: " << std::endl;
	VecDoub a1(10);
	VecDoub b1(10);
	VecDoub c1(10);
	VecDoub d1(10);
	a1.assign(10, 0);
	b1.assign(10, 0);
	c1.assign(10, 0);
	d1.assign(10, 0);

	VecDoub N1(10);
	for (int i = 0; i < 10; ++i) {
		N1[i] = pow(2, i) + 1;
		a1[i] = NewtonsCotes::trapezoidalMethod(f1, N1[i], 0, 1);
		b1[i] = NewtonsCotes::trapezoidalMethod(f2, N1[i], 0, 1);
		c1[i] = NewtonsCotes::trapezoidalMethod(f3, N1[i], 0, 1);
		d1[i] = NewtonsCotes::trapezoidalMethod(f4, N1[i], 0, 1);
	}
	std::cout << "f1: " << std::endl;
	NewtonsCotes::printoutTable(a1, N1);
	std::cout << "f2: " << std::endl;
	NewtonsCotes::printoutTable(b1, N1);
	std::cout << "f3: " << std::endl;
	NewtonsCotes::printoutTable(c1, N1);
	std::cout << "f4: " << std::endl;
	NewtonsCotes::printoutTable(d1, N1);

	std::cout << "Simpson Method: " << std::endl;
	VecDoub a2(10);
	VecDoub b2(10);
	VecDoub c2(10);
	VecDoub d2(10);
	a2.assign(10, 0);
	b2.assign(10, 0);
	c2.assign(10, 0);
	d2.assign(10, 0);
	VecDoub N2(10);
	for (int i = 0; i < 10; ++i) {
		N2[i] = pow(2, i) + 1;
		a2[i] = NewtonsCotes::simpsonsMethod(f1, N2[i], 0, 1);
		b2[i] = NewtonsCotes::simpsonsMethod(f2, N2[i], 0, 1);
		c2[i] = NewtonsCotes::simpsonsMethod(f3, N2[i], 0, 1);
		d2[i] = NewtonsCotes::simpsonsMethod(f4, N2[i], 0, 1);
	}
	std::cout << "f1: " << std::endl;
	NewtonsCotes::printoutTable(a2, N2);
	std::cout << "f2: " << std::endl;
	NewtonsCotes::printoutTable(b2, N2);
	std::cout << "f3: " << std::endl;
	NewtonsCotes::printoutTable(c2, N2);
	std::cout << "f4: " << std::endl;
	NewtonsCotes::printoutTable(d2, N2);

	return 0;
}
*/

//Lektion 9       -----------Ved ikke lige om jeg fik lavet richardson Extrapolation ------------------------
/*
#include <iostream>
#include <cmath>
#include "..\..\NR_C301\code\nr3.h" 
#include "..\..\NR_C301\code\quadrature.h"
#include "..\..\NR_C301\code\derule.h"

using namespace std;

double f1(double x, double tau) {
	double ans = cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f2(double x, double tau) {
	double ans = sqrt(x) * cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f3(double x, double tau) {
	double ans = (1 / sqrt(x)) * cos(pow(x, 2)) * exp(-x);
	return ans;
}

double f4(double x, double tau) {
	double ans = 1000 * exp(-1 / x) * exp(-1 / (1 - x));
	return ans;
}

template<class T>
double derule(T& func, const double a, const double b, const double pre)
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
		std::cout << std::fixed;
		std::cout << std::setprecision(4);
		cout << count << " | " << Ah0 << " | " << abs(Ah0 - Ah1) << endl;
	}
	return Ah0;
}

int main() {
	cout << "DE rule for f1." << endl;
	cout << "Count  |   Ah0  | Ah0-Ah1" << endl;
	derule(f1, 0, 1, 1e-15);
	cout << endl << endl;
	
	cout << "DE rule for f2." << endl;
	cout << "Count  |   Ah0  | Ah0-Ah1" << endl;
	derule(f2, 0, 1, 1e-15);
	cout << endl << endl;

	cout << "DE rule for f3." << endl;
	cout << "Count  |   Ah0  | Ah0-Ah1" << endl;
	derule(f3, 0, 1, 1e-15);
	cout << endl << endl;

	cout << "DE rule for f4." << endl;
	cout << "Count  |   Ah0  | Ah0-Ah1" << endl;
	derule(f4, 0, 1, 1e-15);
	cout << endl << endl;
	
	return 0;
}
*/

//Lektion 10
#include "NumericalSolutions.h"

double Leap_Frog(double x0, double y0, int a, int b, int hinv)
{
	double x = x0;
	double y = y0;

	double t = a;
	double stepsize = (double(b) - double(a)) / double(hinv);

	double xminus1 = x - stepsize * (x * y);		// Using Euler
	double yminus1 = y - stepsize * (-pow(x, 2));	// Using Euler

	//std::cout << "Leap_Frog method, with stepsize " << stepsize << std::endl;
	//std::cout << "x(0) = " << x << std::endl;
	//std::cout << "y(0) = " << y << std::endl;
	//std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
	//std::cout << std::endl;

	for (double i = a; i < b; i = i + stepsize)
	{
		double xplus1 = xminus1 + 2 * stepsize * (x * y);
		double yplus1 = yminus1 + 2 * stepsize * (-pow(x, 2));
		xminus1 = x;
		yminus1 = y;
		x = xplus1;
		y = yplus1;
		t = t + stepsize;
		//	std::cout << "x(" << t << ") = " << x << std::endl;
		//	std::cout << "y(" << t << ") = " << y << std::endl;
		//	std::cout << "x^2 + y^2 = " << pow(x,2)+pow(y,2) << std::endl;
		//	std::cout << std::endl;

	}
	return(x);

}

VecDoub f(double X0, double Y0) {
	VecDoub f(2);
	f[0] = X0 * Y0;
	f[1] = (-pow(X0, 2));
	return f;
}

int main() {
	VecDoub eul(10);
	VecDoub LF(10);
	VecDoub Mp(10);
	VecDoub h(10);
	for (int i = 0; i < eul.size(); ++i) {
		h[i] = 5 * pow(2, i);
	}
	std::cout << "Eulers method: " << std::endl;
	NumericalSolutions::Euler(f, 1, 1, 0, 5, h, eul);
	NumericalSolutions::printoutTable(eul, h);

	std::cout << "Leap Frog method: " << std::endl;
	NumericalSolutions::LeapFrog(f, 1, 1, 0, 5, h, LF);
	NumericalSolutions::printoutTable(LF, h);

	std::cout << "Midpoint method: " << std::endl;
	NumericalSolutions::Midpoint(f, 1, 1, 0, 5, h, Mp);
	NumericalSolutions::printoutTable(Mp, h);



	Doub s, olds = 0.0, oldolds = 0.0, rn = 0.0, r = 0.0, oldr = 0.0, oldoldr = 0.0;
	cout << setw(15) << "hinv" << setw(15) << "Sn" << setw(15) << "k1" << setw(15) << "Rn" << setw(15) << "k2" << setw(15) << "R2n" << endl;
	for (int j = 0; j < 15; j++) {


		//s = Euler(1,1,0,5,5*pow(2,j));  // 1(1) and 2(3)
		s = Leap_Frog(1, 1, 0, 10, 5 * pow(2, j));//1(1) and 2(3)

		//s =Midpoint(1,1,0,5,5*pow(2,j)); // 2(3) and 3(8)

		//s=Trapezoidal(1,1,0,5,5*pow(2,j));// 2(3) and 3(8)

		//s = CalcRK4(1,1,0,5,5*pow(2,j)); //4(15) and 5(24)

		cout << setw(15) << 5 * pow(2, j);
		cout << setw(15) << s;
		if (oldolds != 0.0) { cout << setw(15) << log2(abs((oldolds - olds) / (olds - s))); }
		else { cout << setw(15) << " "; }
		if (olds != 0.0) { r = s + (s - olds) / 1.0; cout << setw(15) << r; } // We see that k1 tends towards 2,and with alpha=2 we get alpha^k-1 = 3
		else { cout << setw(15) << " "; }
		if (oldoldr != 0.0) { cout << setw(15) << log2(abs((oldoldr - oldr) / (oldr - r))); }
		else { cout << setw(15) << " "; }
		if (oldr != 0.0) { rn = r + (r - oldr) / 3.0; cout << setw(15) << rn; } // We see that k2 tends towards 4,and with alpha=2 we get alpha^k-1 = 15
		else { cout << setw(15) << " "; }
		cout << endl;


		oldoldr = oldr;
		oldolds = olds;
		oldr = r;
		olds = s;


	}

	return 0;
}


//Lektion 11
/*

int main() {


	return 0;
}
*/

//Lektion 12
/*
#include <iostream>
#include <cmath>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\nr3.h> 
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\utilities.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\utilities.h>

double F(double yp, double y, double x) {
	return 2 * x + sin(yp) - cos(y);
}

double Fy(double yp, double y, double x) {
	return sin(y);
}

double Fyp(double yp, double y, double x) {
	return cos(yp);
}

std::vector<VecDoub> tridagDiagonals(double N, double h, int a, int b, VecDoub& yi, VecDoub& xi) {
	//Initializer the different vectors for tridag
	VecDoub tA(N-1);
	VecDoub tB(N-1);
	VecDoub tC(N-1);

	tA[0] = 0;
	tB[0] = 2 + pow(h, 2) * Fy(((yi[2] - a) / (2 * h)), yi[1], xi[1]);
	tC[0] = -1 + (h / 2) * Fyp(((yi[2] - a) / (2 * h)), yi[1], xi[1]);
	tA[N - 1] = -1 - (h / 2) * Fyp(((b - yi[N - 2]) / (2 * h)), yi[N - 1], xi[N - 1]);
	tB[N - 1] = 2 + pow(h, 2) * Fy(((b - yi[N - 2]) / (2 * h)), yi[N - 1], xi[N - 1]);

	for (int i = 1; i < N - 2; ++i) {
		tA[i] = -1 - (h / 2) * Fyp(((yi[i + 1] - yi[i - 1]) / (2 * h)), yi[i], xi[i]);
		tB[i] = 2 + pow(h, 2) * Fy(((yi[i + 1] - yi[i - 1]) / (2 * h)), yi[i], xi[i]);
		tC[i] = -1 - (h / 2) * Fyp(((yi[i + 1] - yi[i - 1]) / (2 * h)), yi[i], xi[i]);
	}

	//Creates output vector
	std::vector<VecDoub> finished(3);
	finished.push_back(tA);
	finished.push_back(tB);
	finished.push_back(tC);
	return finished;
}

int main() {
	int a = 0;	//Lower bound
	int b = 2;	//Upper bound
	double N = 3;	//Steps
	double h = (b - a) / (N - 1); //stepsize
	//Intansiere initial guess til en lige linje mellem de to punkter
	VecDoub yi(N);
	VecDoub xi(N);
	yi[0] = a;
	xi[0] = a;
	for (int i = 1; i < N; ++i) {
		yi[i] = yi[i-1] + h;
	}
	for (int i = 1; i < N; ++i) {
		xi[i] = a + i * h;
	}

	std::vector<VecDoub> test(3);
	test = tridagDiagonals(N, h, a, b, yi, xi);
	util::print(test[1]);


	return 0;
}*/

