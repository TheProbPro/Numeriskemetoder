#ifndef NUMERICALSOLUTIONS_H
#define NUMERICALSOLUTIONS_H

// Includes
#include <iostream>
#include <cmath>
#include "../../NR_C301/code/nr3.h"


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
	
	// Eulers method
	VecDoub Euler(VecDoub (*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub& h, VecDoub& ans) {
		VecDoub Fxy(2);
		VecDoub FxyNew(2);
		VecDoub xy(2);
		double stepsize = 0;

		
		for (int i = 0; i < ans.size(); ++i) {
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);
			for (double j = a; j < b; j = j + stepsize) {
				Fxy = f(xy[0], xy[1]);
				FxyNew[0] = xy[0] + stepsize * Fxy[0];
				FxyNew[1] = xy[1] + stepsize * Fxy[1];
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
				ans[i] = xy[0];
			}
		}
		return ans;
	}
	
	//--------------------------------Ikke færdige--------------------------------------------
	// Leap Frog Method
	VecDoub LeapFrog(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub& h, VecDoub& ans) {
		VecDoub Fxy(2);
		VecDoub FxyNew(2);
		VecDoub FxyOld(2);
		VecDoub xy(2);
		double stepsize = 0;

		for (int i = 0; i < ans.size(); ++i) {
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);
			Fxy = f(xy[0], xy[1]);
			FxyOld[0] = xy[0] - stepsize * Fxy[0];
			FxyOld[1] = xy[1] - stepsize * Fxy[1];
			for (double j = a; j < b; j = j + stepsize) {
				Fxy = f(xy[0], xy[1]);
				FxyNew[0] = FxyOld[0] + 2 * stepsize * Fxy[0];
				FxyNew[1] = FxyOld[1] + 2 * stepsize * Fxy[1];
				FxyOld[0] = xy[0];
				FxyOld[1] = xy[1];
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
				ans[i] = xy[0];
			}
		}
		return ans;
	}

	
	// Midpoint method
	VecDoub Midpoint(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		VecDoub Fxy0(2);
		VecDoub Fxy1(2);
		VecDoub FxyNew(2);
		VecDoub FxyHalf(2);
		VecDoub xy(2);
		double stepsize = 0;
		
		for (int i = 0; i < ans.size(); ++i) {
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);
			for (double j = a; j < b; j = j + stepsize) {
				Fxy0 = f(xy[0], xy[1]);
				Fxy1 = f(FxyHalf[0], FxyHalf[1]);
				FxyHalf[0] = xy[0] + 0.5 * stepsize * Fxy0[0];
				FxyHalf[1] = xy[1] + 0.5 * stepsize * Fxy0[1];
				FxyNew[0] = xy[0] + stepsize * Fxy1[0];
				FxyNew[1] = xy[1] + stepsize * Fxy1[1];
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
				ans[i] = xy[0];
			}
		}

		return ans;
	}

	// Trapezoidal method
	VecDoub Trapezoidal(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		VecDoub Fxy(2);
		VecDoub FxyEuler(2);
		VecDoub FxyNew(2);
		VecDoub xy(2);
		double stepsize = 0;

		for (int i = 0; i < ans.size(); ++i) {
			xy[0] = X0;
			xy[1] = Y0;
			stepsize = (double(b) - double(a)) / double(h[i]);
			for (double j = a; j < b; j = j + stepsize) {
				Fxy = f(xy[0], xy[1]);
				
				
				xy[0] = FxyNew[0];
				xy[1] = FxyNew[1];
				ans[i] = xy[0];
			}
		}

		return ans;
	}

	// 4. order Runge Kutta method
	VecDoub RungeKutta(VecDoub(*f)(double X0, double Y0), double X0, double Y0, int a, int b, VecDoub h, VecDoub& ans) {
		return ans;
	}
	
};

#endif // NUMERICALSOLUTIONS_H