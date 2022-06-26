#ifndef NUMERICALSOLUTIONS_H
#define NUMERICALSOLUTIONS_H

// Includes
#include <iostream>
#include <cmath>
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/rk4.h"

namespace NumericalSolutions
{
	//Print method
	void printoutTable(VecDoub& a, VecDoub& fcalc) {

		std::cout << std::left;
		std::cout << std::fixed;
		std::cout << std::setprecision(7);

		std::cout << std::setw(8) << "i" << " | " << std::setw(8) << "a[i]" << " | " <<
			std::setw(8) << "a[i-1] - a[i]" << " | " <<
			std::setw(8) << "Rich-alp^k" << " | " <<
			std::setw(8) << "Rich-fejl" << " | " <<
			std::setw(8) << "Antal f-beregninger" << std::endl << std::endl;

		for (size_t i = 0; i < a.size(); i++) {
			if (i == 0) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else if (i == 1) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << a[i] - a[i - 1] << "| " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << a[i] - a[i - 1] << "  |  " <<
					std::setw(8) << abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) << " | " <<
					std::setw(8) << abs((a[i] - a[i - 1]) / ((abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) - 1))) << " | " <<
					std::setw(8) << fcalc[i] << std::endl << std::endl;
			}
		}

	}

	// Print out method to a given accuracy
	void printoutTable(VecDoub& a, VecDoub& fcalc, double accuracy) {

		std::cout << std::left;
		std::cout << std::fixed;
		std::cout << std::setprecision(7);

		std::cout << std::setw(8) << "i" << " | " << std::setw(8) << "a[i]" << " | " <<
			std::setw(8) << "a[i-1] - a[i]" << " | " <<
			std::setw(8) << "Rich-alp^k" << " | " <<
			std::setw(8) << "Rich-fejl" << " | " <<
			std::setw(8) << "Antal f-beregninger" << std::endl << std::endl;

		for (size_t i = 0; i < a.size(); i++) {
			if (i == 0) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else if (i == 1) {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << a[i] - a[i - 1] << "| " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << " ******* " << " | " <<
					std::setw(8) << fcalc[i] << std::endl;
			}
			else if (i > 6 && (a[i] - a[i - 1]) < accuracy){
				std::cout << "error: " << a[i] - a[i-1] << std::endl;
				break;
			}
			else {
				std::cout << std::setw(8) << i << " | " <<
					std::setw(8) << a[i] << " | " <<
					std::setw(8) << a[i] - a[i - 1] << "  |  " <<
					std::setw(8) << abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) << " | " <<
					std::setw(8) << abs((a[i] - a[i - 1]) / ((abs(a[i - 2] - a[i - 1]) / (abs(a[i - 1] - a[i])) - 1))) << " | " <<
					std::setw(8) << fcalc[i] << std::endl << std::endl;
			}
		}

	}
	
	// Eulers method
	VecDoub Euler(VecDoub (*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub& h, VecDoub& ans) {
		// Declaration of variables
		VecDoub Fxy(2), FxyNew(2), xy(2);
		double stepsize = 0;
		
		for (int i = 0; i < ans.size(); ++i) {
			// Set x and y start values and stepsize for each h value
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);

			for (double j = a; j < b; j = j + stepsize) {
				// Calculate the function for each iteration of x and y values
				Fxy = f(xy[0], xy[1]);
				FxyNew[0] = xy[0] + stepsize * Fxy[0];
				FxyNew[1] = xy[1] + stepsize * Fxy[1];
				// Setting the new x and y values
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
			}
			// Saving the x values to the answer vector
			ans[i] = xy[0];
		}
		return ans;
	}
	
	// Leap Frog Method
	VecDoub LeapFrog(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub& h, VecDoub& ans) {
		// Declaration of variables
		VecDoub Fxy(2), FxyNew(2), FxyOld(2), xy(2);
		double stepsize = 0;

		for (int i = 0; i < ans.size(); ++i) {
			// Set x and y start values and stepsize for each h value
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);
			// Calculate yn-1 using euler
			Fxy = f(xy[0], xy[1]);
			FxyOld[0] = xy[0] - stepsize * Fxy[0];
			FxyOld[1] = xy[1] - stepsize * Fxy[1];

			for (double j = a; j < b; j = j + stepsize) {
				// Calculate the function for each iteration of x and y values
				Fxy = f(xy[0], xy[1]);
				FxyNew[0] = FxyOld[0] + 2 * stepsize * Fxy[0];
				FxyNew[1] = FxyOld[1] + 2 * stepsize * Fxy[1];
				// Setting the new x and y values and updating the old ones
				FxyOld[0] = xy[0];
				FxyOld[1] = xy[1];
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
			}
			// Saving the x values to the answer vector
			ans[i] = xy[0];
		}
		return ans;
	}

	
	// Midpoint method
	VecDoub Midpoint(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		// Declaration of variables
		VecDoub Fxy0(2), Fxy1(2), FxyNew(2), FxyHalf(2), xy(2);
		double stepsize = 0;
		
		for (int i = 0; i < ans.size(); ++i) {
			// Set x and y start values and stepsize for each h value
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);

			for (double j = a; j < b; j = j + stepsize) {
				// Calculate the function for each iteration of x and y values
				Fxy0 = f(xy[0], xy[1]);
				FxyHalf[0] = xy[0] + 0.5 * stepsize * Fxy0[0];
				FxyHalf[1] = xy[1] + 0.5 * stepsize * Fxy0[1];
				Fxy1 = f(FxyHalf[0], FxyHalf[1]);
				FxyNew[0] = xy[0] + stepsize * Fxy1[0];
				FxyNew[1] = xy[1] + stepsize * Fxy1[1];
				// Setting the new x and y values
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
			}
			// Saving the x values to the answer vector
			ans[i] = xy[0];
		}
		return ans;
	}

	// Trapezoidal method
	VecDoub Trapezoidal(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		// Declaration of variables
		VecDoub Fxy0(2), Fxy1(2), FxyEuler(2), FxyNew(2), xy(2);
		double stepsize = 0;

		for (int i = 0; i < ans.size(); ++i) {
			// Set x and y start values and stepsize for each h value
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);

			for (double j = a; j < b; j = j + stepsize) {
				// Calculate the function for each iteration of x and y values
				Fxy0 = f(xy[0], xy[1]);
				FxyEuler[0] = xy[0] + stepsize * Fxy0[0];
				FxyEuler[1] = xy[1] + stepsize * Fxy0[1];

				Fxy1 = f(FxyEuler[0], FxyEuler[1]);
				FxyNew[0] = xy[0] + 0.5 * stepsize * (Fxy1[0] + Fxy0[0]);
				FxyNew[1] = xy[1] + 0.5 * stepsize * (Fxy1[1] + Fxy0[1]);
				// Setting the new x and y values
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
			}
			// Saving the x values to the answer vector
			ans[i] = xy[0];
		}
		return ans;
	}

	// 4. order Runge Kutta method
	VecDoub RungeKutta(VecDoub(*f)(double X0, double Y0), void(*dydx)(const Doub x, VecDoub_I& y, VecDoub_O& dydx), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		// Set x and y start values and stepsize for each h value
		VecDoub xy(2), dydx1(2), xyOld(2);
		double stepsize = 0;

		for (int i = 0; i < ans.size(); ++i) {
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);

			for (double j = a; j < b; j = j + stepsize) {
				dydx(a, xy, dydx1);
				xyOld = xy;
				rk4(xyOld, dydx1, a, stepsize, xy, dydx);
			}
			ans[i] = xy[0];
		}
		return ans;
	}
	
};

#endif // NUMERICALSOLUTIONS_H