#ifndef NEWTONSMETHOD_H
#define NEWTONSMETHOD_H

#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/qrdcmp.h"
#include "../../NR_C301/code/roots_multidim.h"

namespace newtonsMethod {
	template <class T>
	void newt(VecDoub_IO& x, Bool& check, T& vecfunc, double iterations) {
		const Int MAXITS = iterations;
		const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
		const Doub TOLX = numeric_limits<Doub>::epsilon();
		Int i, j, its, n = x.size();
		Doub den, f, fold, stpmax, sum, temp, test;
		VecDoub g(n), p(n), xold(n);
		MatDoub fjac(n, n);
		NRfmin<T> fmin(vecfunc);
		NRfdjac<T> fdjac(vecfunc);
		VecDoub& fvec = fmin.fvec;
		f = fmin(x);
		test = 0.0;
		for (i = 0; i < n; i++)
			if (abs(fvec[i]) > test) test = abs(fvec[i]);
		if (test < 0.01 * TOLF) {
			check = false;
			return;
		}
		sum = 0.0;
		for (i = 0; i < n; i++) sum += SQR(x[i]);
		stpmax = STPMX * MAX(sqrt(sum), Doub(n));
		for (its = 0; its < MAXITS; its++) {
			fjac = fdjac(x, fvec);
			for (i = 0; i < n; i++) {
				sum = 0.0;
				for (j = 0; j < n; j++) sum += fjac[j][i] * fvec[j];
				g[i] = sum;
			}
			for (i = 0; i < n; i++) xold[i] = x[i];
			fold = f;
			for (i = 0; i < n; i++) p[i] = -fvec[i];
			LUdcmp alu(fjac);
			alu.solve(p, p);
			lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
			test = 0.0;
			for (i = 0; i < n; i++)
				if (abs(fvec[i]) > test) test = abs(fvec[i]);
			if (test < TOLF) {
				check = false;
				return;
			}
			if (check) {
				test = 0.0;
				den = MAX(f, 0.5 * n);
				for (i = 0; i < n; i++) {
					temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
					if (temp > test) test = temp;
				}
				check = (test < TOLMIN);
				return;
			}
			test = 0.0;
			for (i = 0; i < n; i++) {
				temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
				if (temp > test) test = temp;
			}
			if (test < TOLX)
				return;
		}
		std::cout << "MAXITS exceeded in newt" << std::endl;
	}

	template <class T>
	void newt(VecDoub_IO& x, Bool& check, T& vecfunc, double iterations, double accuracy) {
		const Int MAXITS = iterations;
		const Doub TOLF = accuracy, TOLMIN = 1.0e-12, STPMX = 100.0;
		const Doub TOLX = numeric_limits<Doub>::epsilon();
		Int i, j, its, n = x.size();
		Doub den, f, fold, stpmax, sum, temp, test;
		VecDoub g(n), p(n), xold(n);
		MatDoub fjac(n, n);
		NRfmin<T> fmin(vecfunc);
		NRfdjac<T> fdjac(vecfunc);
		VecDoub& fvec = fmin.fvec;
		f = fmin(x);
		test = 0.0;
		for (i = 0; i < n; i++)
			if (abs(fvec[i]) > test) test = abs(fvec[i]);
		if (test < 0.01 * TOLF) {
			check = false;
			return;
		}
		sum = 0.0;
		for (i = 0; i < n; i++) sum += SQR(x[i]);
		stpmax = STPMX * MAX(sqrt(sum), Doub(n));
		for (its = 0; its < MAXITS; its++) {
			fjac = fdjac(x, fvec);
			for (i = 0; i < n; i++) {
				sum = 0.0;
				for (j = 0; j < n; j++) sum += fjac[j][i] * fvec[j];
				g[i] = sum;
			}
			for (i = 0; i < n; i++) xold[i] = x[i];
			fold = f;
			for (i = 0; i < n; i++) p[i] = -fvec[i];
			LUdcmp alu(fjac);
			alu.solve(p, p);
			lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
			test = 0.0;
			for (i = 0; i < n; i++)
				if (abs(fvec[i]) > test) test = abs(fvec[i]);
			if (test < TOLF) {
				check = false;
				return;
			}
			if (check) {
				test = 0.0;
				den = MAX(f, 0.5 * n);
				for (i = 0; i < n; i++) {
					temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
					if (temp > test) test = temp;
				}
				check = (test < TOLMIN);
				return;
			}
			test = 0.0;
			for (i = 0; i < n; i++) {
				temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
				if (temp > test) test = temp;
			}
			if (test < TOLX)
				return;
		}
		std::cout << "MAXITS exceeded in newt" << std::endl;
	}
}


#endif // NEWTONSMETHOD