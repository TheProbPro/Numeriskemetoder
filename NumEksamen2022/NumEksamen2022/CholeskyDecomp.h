#ifndef CHOLESKYDECOMP_H
#define CHOLESKYDECOMP_H

/*
	This code is inspired by code i made during the lectures.
	During the lectures i have talked to Simon Vinkel, Thomas Therkelsen, Mark Egedal, Anders Lind-Thomsen og Peter Frydensberg
*/

//Includes
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/cholesky.h"
#include "../../NR_C301/code/utilities.h"

namespace CholeskyDecomp
{
	void CholeskyDecomp(MatDoub& C, VecDoub& c, int size) {
		// Createing vector to contain answer
		VecDoub chCoef(c.size());

		// Solving cholesky
		Cholesky ch(C);
		ch.solve(c, chCoef);

		// Finding L and pinting the diagonal
		std::cout << "Diagonal elements of L: " << std::endl;
		MatDoub L = ch.el;
		VecDoub LDiag(6, 0.0);
		for (int i = 0; i < L.ncols(); ++i) {
			LDiag[i] = L[i][i];
		}
		util::print(LDiag);
		std::cout << std::endl;

		// Printing answer
		std::cout << "Cholesky Decomposition answers: " << std::endl;
		util::print(chCoef);
	}
};

#endif // CHOLESKYDECOMP_H