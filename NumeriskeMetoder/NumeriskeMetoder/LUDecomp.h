#ifndef LUDECOMP_H
#define LUDECOMP_H

//Includes
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/utilities.h"

namespace LUDecomp {
	void LUDecomp(MatDoub& A, VecDoub& b, int size) {
		// Evaluate x
		VecDoub x(size);
		//Create LU decomposition object
		LUdcmp LU(A);

		//Create LU matrix and print
		auto L = LU.lu;
		for (int i = 0; i < size-1; ++i) {
			for (int j = 1; j < size; ++j) {
				if (i == j) {
					;
				}
				else {
					L[i][j] = 0;
				}
			}
		}
		for (int i = 0; i < size; ++i) {
			L[i][i] = 1;
		}
		util::print(L, "L");

		auto U = LU.lu;
		for (int i = 1; i < size; ++i) {
			for (int j = 0; j < size-1; ++j) {
				if (i == j) {
					;
				}
				else {
					U[i][j] = 0;
				}
			}
		}
		util::print(U, "U");

		LU.solve(b, x);

		// print x
		util::print(A, "A");
		util::print(L * U, "L*U");
		util::print(L * U * x, "L*U*x");
		util::print(x, "x");
	}

	void LUDecompquick(MatDoub& A, VecDoub& b, int size) {
		LUdcmp lu(A);
		VecDoub ans(size);

		lu.solve(b, ans);
		std::cout << "LU Decomposition answer: " << std::endl;
		util::print(ans);
	}
};


#endif // LUDECOMP_H