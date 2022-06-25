#ifndef CONVERGENCEMETHODS_H
#define CONVERGENCEMETHODS_H

// Includes
#define _USE_MATH_DEFINES
#include <math.h>
#include "..\..\NR_C301\code\nr3.h"
#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/qrdcmp.h"
#include "../../NR_C301/code/roots_multidim.h"
#include "..\..\NR_C301\code\utilities.h"
#include <vector>

namespace ConvergenceMethods
{
	//Print metode
	void printOutput(vector<double> x, double k /* order */) {
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

	//Newtons method
	void Newton(double (*f)(double x), double (*df)(double x), int iterations) {
		std::vector<double> x;
		x.push_back(0);
		for (int i = 0; i < iterations; ++i) {
			x.push_back(x[i] - f(x[i]) / df(x[i]));
		}
		printOutput(x, 2);

		cout << "Estimated error: " << abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 2)) * pow(x[iterations] - x[iterations - 1], 2) << endl << endl;
	}

	//Secant method

	void Secant(double (*f)(double x), double (*df)(double x), int iterations) {
		vector<double> x;
		x.push_back(0);
		x.push_back(M_PI_2);
		for (int i = 1; i < iterations; ++i) {
			x.push_back(x[i] - ((x[i] - x[i - 1]) / (f(x[i]) - f(x[i - 1])) * f(x[i])));
		}
		printOutput(x, 1.62);

		cout << "Estimated error: " << abs(x[iterations] - x[iterations - 1]) / (pow(abs(x[iterations - 1] - x[iterations - 2]), 1.62)) * pow(x[iterations] - x[iterations - 1], 1.62) << endl << endl;
	}

	// Bitsection method
	void Bitsection(double (*f)(double x), double (*df)(double x), int iterations) {
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

	// Regula Falsi method
	void RegulaFalsi(double (*f)(double x), double (*df)(double x), int iterations) {
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

	// Ridders method
	void Ridders(double (*f)(double x), int iterations, double a, double b) {
		vector<double> x, y, z;
		x.push_back(a);
		y.push_back(b);
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
};

#endif // CONVERGENCEMETHODS_H