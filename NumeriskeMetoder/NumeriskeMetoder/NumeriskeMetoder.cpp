
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
#include <iostream>
#include <fstream>
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/ludcmp.h"
# include "../../NR_C301/code/cholesky.h"
#include "../../NR_C301/code/utilities.h"
//Include til lektion 3 opgave
#include "../../NR_C301/code/svd.h"

using namespace std;

int main() {
	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("C:\\Users\\nvigg\\Desktop\\4.semester Boeger\\4.Sem Real\\NumeriskeMetoder\\Programmeringsopgaver\\FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}
	std::cout << "Filip loaded" << std::endl;
	
	VecDoub xPont(40); VecDoub yPont(40);
	ifstream Pont("C:\\Users\\nvigg\\Desktop\\4.semester Boeger\\4.Sem Real\\NumeriskeMetoder\\Programmeringsopgaver\\PontiusData.dat");
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

	LUdcmp lu(CPont);
	VecDoub pontCoef(cPont.size());

	lu.solve(cPont, pontCoef);
	std::cout << "LU Decomposition answer: " << std::endl;
	util::print(pontCoef);

	//Cholesky
	VecDoub pontCoefch(cPont.size());
	Cholesky ch(CPont);
	ch.solve(cPont, pontCoefch);
	std::cout << "Cholesky Decomposition answers: " << std::endl;
	util::print(pontCoefch);


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

	LUdcmp luF(CFilip);
	VecDoub filipCoef(cFilip.size());

	luF.solve(cFilip, filipCoef);
	std::cout << "LU Decomposition answer: " << std::endl;
	util::print(filipCoef);

	//Cholesky
	
	VecDoub filipCoefch(cFilip.size());
	Cholesky chF(CFilip);
	chF.solve(cFilip, filipCoefch);
	std::cout << "Cholesky Decomposition answers: " << std::endl;
	util::print(filipCoefch);
	
	//Dette deviere fra alle de andre svar af en eller anden grund
	
	//Lektion 3 SVD-decomposition af filip og pontius
	SVD svdP(aPont);
	VecDoub xxPont(3);
	svdP.solve(bPont, xxPont, svdP.eps);
	cout << "SVD Pontius: " << endl;
	util::print(xxPont);

	SVD svdF(aFilip);
	VecDoub xxFil(11);
	svdF.solve(bFilip, xxFil, svdF.eps);
	cout << "SVD Filip: " << endl;
	util::print(xxFil);

	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file




//Lektion 5 + Lektion 6
/*
#define _USE_MATH_DEFINES
#include <cmath>
#include "C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\nr3.h"
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\utilities.h>
#include <vector>


using namespace std;

double f(double x) {
	return x - cos(x);
}

double df(double x) {
	return 1 - sin(x);
}


void printOutput(vector<double> x, double k /* order *//*)  {
	cout << left;
	cout << fixed;
	cout << setprecision(7);

	cout << setw(3) << "i" << " | " << setw(8) << "x[i]" << " | " <<
		setw(8) << "x[i] - x[i-1]" << " | " <<
		"|(x[i]-x[i-1])|/|(x[i-1] - x[i-2])|^k" << endl << endl;
	for (size_t i = 0; i < x.size(); i++) {
		if (i == 0) {
			cout << i << " | " << x[i] << endl;
		}
		else if (i == 1) {
			cout << i << " | " << x[i] << " | " <<
				x[i] - x[i - 1] << endl;
		}
		else {
			cout << i << " | " << x[i] << " | " <<
				x[i] - x[i - 1] << "  |  " <<
				abs(x[i] - x[i - 1]) / (pow(abs(x[i - 1] - x[i - 2]), k)) << endl << endl;
		}
	}
}


void Newton(int iterations) {
	std::vector<double> x;
	x.push_back(0);
	for (int i = 0; i < iterations; ++i) {
		x.push_back(x[i] - f(x[i]) / df(x[i]));
	}
	printOutput(x, 2);

	cout << "Estimated error: " << abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 2)) * pow(x[iterations] - x[iterations - 1], 2) << endl << endl;
}

void Secant(int iterations) {
	vector<double> x;
	x.push_back(0);
	x.push_back(M_PI_2);
	for (int i = 1; i < iterations; ++i) {
		x.push_back(x[i] - ((x[i] - x[i - 1]) / (f(x[i]) - f(x[i - 1])) * f(x[i])));
	}
	printOutput(x, 1.62);

	cout << "Estimated error: " << abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 1.62)) * pow(x[iterations] - x[iterations - 1], 1.62) << endl << endl;
}

void Bitsection(int iterations) {
	vector<double> x, y;
	x.push_back(0);
	y.push_back(M_PI_2);
	for (int i = 0; i < iterations; ++i) {
		x.push_back(((x[i] + y[i]) / 2));
		if (f(x[i + 1]) * f(y[i]) < 0) {
			y.push_back(y[i]);
		}
		else {
			y.push_back(x[i]);
		}
	}
	printOutput(x, 1);

	cout << "Estimated error: " << x[iterations] - x[iterations - 1] << endl << endl;
}

void RegulaFalsi(int iterations) {
	vector<double> x, y;
	x.push_back(0);
	y.push_back(M_PI_2);
	for (int i = 0; i < iterations; ++i) {
		x.push_back(x[i] - ((x[i] - y[i]) / (f(x[i]) - f(y[i]))) * f(x[i]));
		if (f(x[i + 1]) * f(y[i]) < 0) {
			y.push_back(y[i]);
		}
		else {
			y.push_back(x[i]);
		}
	}
	printOutput(x, 1);
	
	double c = abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 1));
	cout << "Estimated error: " << ((-c) / (1 - c)) * (x[iterations] - x[iterations - 1]) << endl << endl;
}

void Ridders(int iterations) {
	vector<double> x, y, z;
	x.push_back(0);
	y.push_back(M_PI_2);
	for (int i = 0; i < iterations; ++i) {
		z.push_back((x[i] + y[i]) / 2);
		bool sign = false;
		if (f(x[i]) > f(y[i])) {
			sign = true;
		}
		if (sign == true) {
			x.push_back(z[i] + (z[i] - x[i]) * ((f(z[i])) / (sqrt(pow(f(z[i]), 2) - (f(x[i]) * f(y[i]))))));
		}
		else if (sign == false) {

			x.push_back(z[i] + (z[i] - x[i]) * ((-f(z[i])) / (sqrt(pow(f(z[i]), 2) - (f(x[i]) * f(y[i]))))));
		}
		if (f(x[i + 1]) * f(z[i]) < 0) {
			y.push_back(z[i]);
		}
		else if (f(x[i + 1]) * f(y[i]) < 0) {
			y.push_back(y[i]);
		}
		else {
			y.push_back(x[i]);
		}
	}
	printOutput(x, 3);
	
	cout << "Estimated error: " << abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 3)) * pow(x[iterations] - x[iterations - 1], 3) << endl << endl;
}



int main() {

	//Newtons method
	cout << "Newtons Method:" << endl;
	Newton(10);
	//Secant method
	cout << "Secant Method:" << endl;
	Secant(8);
	//Bitsection method
	cout << "Bitsection Method:" << endl;
	Bitsection(10);
	//False position method
	cout << "False Position (Regula Falsi) Method:" << endl;
	RegulaFalsi(10);
	//Ridders method
	cout << "Ridders Method:" << endl;
	Ridders(5);


	return 0;
}
*/


/*
//Lektion 7
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\nr3.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\ludcmp.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\qrdcmp.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\roots_multidim.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\utilities.h>

double n = 5;

void printOutput(MatDoub& x, double n) {
	//x[0] = L0, x[1] = L, x[2] = p, x[3] = x, x[4] = thetha, x[5] = phi, x[6] = a, x[7] = H
	cout << left;
	cout << fixed;
	cout << setprecision(7);

	cout << setw(8) << "n" << " | " << setw(12) << "L0" << " | " << setw(12) << "L" << " | " <<
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

//Lektion 8
/*
#include <iostream>
#include <cmath>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\nr3.h> 


void printoutTable(VecDoub& a, VecDoub& fcalc) {
	
	cout << left;
	cout << fixed;
	cout << setprecision(7);

	cout << setw(8) << "i" << " | " << setw(8) << "a[i]" << " | " <<
		setw(8) << "a[i-1] - a[i]" << " | " << 
		setw(8) << "Rich-alp^k" << " | " << 
		setw(8) << "Antal f-beregninger" << endl << endl;
	
	for (size_t i = 0; i < a.size(); i++) {
		if (i == 0) {
			cout << setw(8) << i << " | " <<
				setw(8) << a[i] << " | " << 
				setw(8) << " *** " << " | " <<
				setw(8) << " *** " << " | " <<
				setw(8) << fcalc[i] << endl;
		}
		else if (i == 1) {
			cout << setw(8) << i << " | " << 
				setw(8) << a[i] << " | " <<
				setw(8) << a[i] - a[i - 1] << "| " <<
				setw(8) << " *** " << " | " <<
				setw(8) << fcalc[i] << endl;
		}
		else {
			cout << setw(8) << i << " | " << 
				setw(8) << a[i] << " | " <<
				setw(8) << a[i] - a[i - 1] << "  |  " <<
				setw(8) << abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) << " | " <<
				setw(8) << fcalc[i] << endl << endl;
		}
	}

}

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

double rectangleMethod(double (*f)(double x), double N /*Nummer af inddelinger*//*, double a /*start værdi*//*, double b /*slut værdi*//*) {
	double h = (b - a) / (N - 1); //Stepsize
	double sum = 0;
	for (int i = 0; i < N - 1; ++i) {
		sum += f(a + h * i + h * 0.5);
	}
	double result = h * sum;
	return result;
}

double trapezoidalMethod(double (*f)(double x), double N /*Nummer af inddelinger*//*, double a /*start værdi*//*, double b /*slut værdi*//*) {
	double h = (b - a) / (N - 1); //Stepsize
	double sum = 1/2*f(a)+1/2*f(b);
	for (int i = 1; i < N - 2; ++i) {
		sum += f(a + h * i);
	}
	double result = h * sum;
	return result;
}

double simpsonsMethod(double (*f)(double x), double N /*Nummer af inddelinger*//*, double a /*start værdi*//*, double b /*slut værdi*//*) {
	double h = (b - a) / (N - 1); //Stepsize
	double sum = 1 / 3 * f(a) + 1 / 3 * f(b);
	for (int i = 1; i < N - 2; i+2) {
		sum += 4 / 3 * f(a + h * i);
	}
	for (int i = 2; i < N - 2; i + 2) {
		sum += 2 / 3 * f(a + h * i);
	}
	double result = h * sum;
	return result;
}

int main() {
	VecDoub a(10);
	VecDoub b(10);
	VecDoub c(10);
	VecDoub d(10);

	VecDoub N(10);
	for (int i = 0; i < 10; ++i) {
		N[i] = pow(2, i) + 1;
		a[i] = rectangleMethod(f1, N[i], 0, 1);
		b[i] = rectangleMethod(f2, N[i], 0, 1);
		c[i] = rectangleMethod(f3, N[i], 0, 1);
		d[i] = rectangleMethod(f4, N[i], 0, 1);
	}
	std::cout << "Rectangle Method: " << std::endl;
	std::cout << "f1: " << std::endl;
	printoutTable(a, N);
	std::cout << "f2: " << std::endl;
	printoutTable(b, N);
	std::cout << "f3: " << std::endl;
	printoutTable(c, N);
	std::cout << "f4: " << std::endl;
	printoutTable(d, N);
	
	std::cout << "Trapezoidal Method: " << std::endl;
	VecDoub a1(10);
	VecDoub b1(10);
	VecDoub c1(10);
	VecDoub d1(10);

	VecDoub N1(10);
	for (int i = 0; i < 10; ++i) {
		N1[i] = pow(2, i) + 1;
		a1[i] = trapezoidalMethod(f1, N1[i], 0, 1);
		b1[i] = trapezoidalMethod(f2, N1[i], 0, 1);
		c1[i] = trapezoidalMethod(f3, N1[i], 0, 1);
		d1[i] = trapezoidalMethod(f4, N1[i], 0, 1);
	}
	std::cout << "f1: " << std::endl;
	printoutTable(a1, N1);
	std::cout << "f2: " << std::endl;
	printoutTable(b1, N1);
	std::cout << "f3: " << std::endl;
	printoutTable(c1, N1);
	std::cout << "f4: " << std::endl;
	printoutTable(d1, N1);

	std::cout << "Simpson Method: " << std::endl;
	VecDoub a2(10);
	VecDoub b2(10);
	VecDoub c2(10);
	VecDoub d2(10);

	VecDoub N2(10);
	for (int i = 0; i < 10; ++i) {
		N2[i] = pow(2, i) + 1;
		a2[i] = rectangleMethod(f1, N2[i], 0, 1);
		b2[i] = rectangleMethod(f2, N2[i], 0, 1);
		c2[i] = rectangleMethod(f3, N2[i], 0, 1);
		d2[i] = rectangleMethod(f4, N2[i], 0, 1);
	}
	std::cout << "Rectangle Method: " << std::endl;
	std::cout << "f1: " << std::endl;
	printoutTable(a2, N2);
	std::cout << "f2: " << std::endl;
	printoutTable(b2, N2);
	std::cout << "f3: " << std::endl;
	printoutTable(c2, N2);
	std::cout << "f4: " << std::endl;
	printoutTable(d2, N2);

	return 0;
}
*/

//Lektion 9       -----------Ved ikke lige om jeg fik lavet richardson Extrapolation ------------------------
/*
#include <iostream>
#include <cmath>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\nr3.h> 
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\quadrature.h>
#include <C:\Users\nvigg\Desktop\4.semester Boeger\4. Semester\NR_C301\code\derule.h>

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
		cout << count << " | " << Ah0 << " | " << abs(Ah0 - Ah1) << endl;
	}
	return Ah0;
}

int main() {
	cout << "DE rule for f1." << endl;
	cout << "count | Ah0 | Ah0-Ah1" << endl;
	derule(f1, 0, 1, 1e-15);
	cout << endl << endl;
	
	cout << "DE rule for f2." << endl;
	cout << "count | Ah0 | Ah0-Ah1" << endl;
	derule(f2, 0, 1, 1e-15);
	cout << endl << endl;

	cout << "DE rule for f3." << endl;
	cout << "count | Ah0 | Ah0-Ah1" << endl;
	derule(f3, 0, 1, 1e-15);
	cout << endl << endl;

	cout << "DE rule for f4." << endl;
	cout << "count | Ah0 | Ah0-Ah1" << endl;
	derule(f4, 0, 1, 1e-15);
	cout << endl << endl;
	
	return 0;
}
*/

//Lektion 10
/*
int main() {
	

	return 0;
}
*/

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
