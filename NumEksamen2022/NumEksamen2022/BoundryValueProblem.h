#ifndef BOUNDRYVALUEPROBLEM_H
#define BOUNDRYVALUEPROBLEM_H

/*
	This code is heavely influenced by the code Ole uploaded to itslearning.
	I have rewrithen it partly to make it a function i could call.
*/

// Includes
#include <iostream>
#include <cmath>
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/tridag.h"
#include "../../NR_C301/code/utilities.h"
//#include "../../NR_C301/code/MatPlot.h"

namespace BoundryValueProblem
{
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

	void solve(double (*F)(Doub x, Doub y, Doub yp), double (*Fy)(Doub x, Doub y, Doub yp), double (*Fyp)(Doub x, Doub y, Doub yp), Doub a, Doub b, Doub aa, Doub bb) {
		int N = 1;
		int kk, kkmax = 10;
		Doub h = (b - a) / N;
		yi yinit(aa, bb, N);
		VecDoub y(N + 1, 0.0);
		for (int i = 0; i < N + 1; i++) { y[i] = yinit(i); } //Startguess for Y, first run
		Doub s, olds = 0.0, oldolds = 0.0, rn = 0.0, r = 0.0, oldr = 0.0, oldoldr = 0.0;
		VecDoub yhalf(kkmax, 0.0);
		// Richardson Loop Starts Here
		for (kk = 0; kk < kkmax; kk++)
		{
			N *= 2;
			Doub h = (b - a) / N;
			VecDoub ytemp(N + 1, 0.0);
			VecDoub xplot(N + 1, 0.0);
			VecDoub dy(N + 1, 0.0);
			VecDoub Jp(N + 1, 0.0); //Note that Jp is one shorter than J, but in Tridag they have the same length and the LAST element is ignored
			VecDoub J(N + 1, 0.0);
			VecDoub Jm(N + 1, 0.0); //Note that Jm is one shorter than J, but in Tridag they have the same length and the FIRST element is ignored
			VecDoub Fi(N + 1, 0.0);
			VecDoub Fm(N + 1, 0.0);
			xi x(a, h);
			for (int i = 0; i < N + 1; i++) {
				if (i % 2 == 0) { ytemp[i] = y[i / 2]; }
				else { ytemp[i] = (y[(i + 1) / 2] + y[(i - 1) / 2]) / 2.0; }
			} //Startguess for Y , second runs
			y = ytemp;
			for (int i = 0; i < N + 1; i++) { xplot[i] = x(i); } //corresponding x values


			//"Newton" LOOP STARTS HERE
			for (int k = 0; k < 3; k++) {
				// Define Fi vector
				Fi[0] = y[0] - aa;
				Fi[N] = y[N] - bb;
				for (int i = 1; i < N; i++) { Fi[i] = y[i + 1] - 2 * y[i] + y[i - 1] - h * h * F(x(i), y[i], (y[i + 1] - y[i - 1]) / (h * 2.0)); }
				for (int i = 0; i < N + 1; i++) { Fm[i] = -Fi[i]; }
				// Define JVectors
				J[0] = 1;
				J[N] = 1;
				for (int i = 1; i < N; i++)
				{
					Jm[i] = 1 + h / 2.0 * Fyp(x(i), y[i], (y[i + 1] - y[i - 1]) / (h * 2.0));        // Define the Ji,i-1 vector, called Jminus or Jm
					J[i] = -2 - h * h * Fy(x(i), y[i], (y[i + 1] - y[i - 1]) / (h * 2.0));        // Define the Ji,i vector, called J
					Jp[i] = 1 - h / 2.0 * Fyp(x(i), y[i], (y[i + 1] - y[i - 1]) / (h * 2.0));        // Define the Ji,i+1 vector, called Jplus or Jp
				}
				// if you want, put them in a Matrix.. alternatively use Tridag


				tridag(Jm, J, Jp, Fm, dy); //Find dy
				for (int i = 0; i < N + 1; i++) { y[i] = y[i] + dy[i]; } // Update y
			}
			//Richardson Calculations
			yhalf[kk] = y[N / 2];
			s = y[N / 2];
			cout << setw(15) << kk + 1;
			cout << setw(15) << s;
			if (abs(oldolds) > 0.0) { cout << setw(15) << log2((oldolds - olds) / (olds - s)); }
			else { cout << setw(15) << " "; }
			if (abs(olds) > 0.0) { r = s + (s - olds) / 3; cout << setw(15) << r; } // We see that k1 tends towards 2,and with alpha=2 we get alpha^k-1 = 3
			else { cout << setw(15) << " "; }
			if (abs(oldoldr) > 0.0) { cout << setw(15) << log2((oldoldr - oldr) / (oldr - r)); }
			else { cout << setw(15) << " "; }
			if (abs(oldr) > 0.0) { rn = r + (r - oldr) / 15; cout << setw(15) << rn; } // We see that k2 tends towards 4,and with alpha=2 we get alpha^k-1 = 15
			else { cout << setw(15) << " "; }
			cout << endl;
			Doub alphak = abs((oldolds - olds) / (olds - s));

			if (kk > 3)
				if (abs(s - olds) / alphak < 0.0001) // Terminate at a predetermined error estimate
					break;

			oldoldr = oldr;
			oldolds = olds;
			oldr = r;
			olds = s;

		}
		VecDoub yprint(kk, 0.0);
		for (int i = 0; i < kk; i++)
		{
			yprint[i] = yhalf[i];
		}
		util::print(yprint, "y at midpoint");
	}
};

#endif // BOUNDRYVALUEPROBLEM_H