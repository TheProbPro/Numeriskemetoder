#ifndef NEWTONSCOTES_H
#define NEWTONSCOTES_H

/*
	This code is inspired by code i made during the lectures.
	During the lectures i have talked to Simon Vinkel, Thomas Therkelsen, Mark Egedal, Anders Lind-Thomsen og Peter Frydensberg
*/

// Includes
#include <iostream>
#include <cmath>
#include "..\..\NR_C301\code\nr3.h"
#include "..\..\NR_C301\code\quadrature.h"
#include "..\..\NR_C301\code\derule.h"

namespace NewtonsCotes
{
	//Print method
	void printoutTable(VecDoub& a, VecDoub& fcalc) {

		std::cout << std::left;
		std::cout << std::fixed;
		std::cout << std::setprecision(5);

		std::cout << std::setw(8) << "i" << " | " << std::setw(8) << "a[i]" << " | " <<
			std::setw(8) << "a[i-1] - a[i]" << " | " <<
			std::setw(8) << "Rich-alp^k" << " | " <<
			std::setw(8) << "Rich-fejl" << " | " <<
			std::setw(8) << "Antal f-beregninger" << std::endl << std::endl;

		for (size_t i = 0; i < a.size(); i++) {
			if (i == 0) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(13) << " ******* " << " | " <<
					std::setw(10) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else if (i == 1) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(14) << a[i] - a[i - 1] << "| " <<
					std::setw(10) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(12) << a[i] - a[i - 1] << "  |  " <<
					std::setw(9) << abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) << " | " <<
					std::setw(9) << abs((a[i] - a[i - 1]) / ((abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) - 1))) << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
		}

	}

	//Rectangle method
	double rectangleMethod(double (*f)(double x), double N /*Nummer af inddelinger*/, double a /*start værdi*/, double b /*slut værdi*/) {
		double h = (b - a) / (N - 1); //Stepsize
		double sum = a;
		for (int i = 0; i < N - 1; ++i) {
			sum += f(a + h * i + h * 0.5);
		}
		double result = h * sum;
		return result;
	}

	//Trapezoidal method
	double trapezoidalMethod(double (*f)(double x), double N /*Nummer af inddelinger*/, double a /*start værdi*/, double b /*slut værdi*/) {
		double h = (b - a) / (N - 1); //Stepsize
		double sum = 0.5 * f(a) + 0.5 * f(b);
		for (int i = 1; i < N - 1; ++i) {
			sum += f(a + h * i);
		}
		double result = h * sum;
		return result;
	}

	//Simpsons method
	double simpsonsMethod(double (*f)(double x), double N /*Nummer af inddelinger*/, double a /*start værdi*/, double b /*slut værdi*/) {
		double h = (b - a) / (N - 1); //Stepsize
		double sum = static_cast<double>(1) / 3 * f(a) + static_cast<double>(1) / 3 * f(b);
		for (int i = 1; i < N - 1; i = i + 2) {
			sum += static_cast<double>(4) / 3 * f(a + h * i);
		}
		for (int i = 2; i < N - 1; i = i + 2) {
			sum += static_cast<double>(2) / 3 * f(a + h * i);
		}
		double result = h * sum;
		return result;
	}

};

#endif // NEWTONSCOTES_H