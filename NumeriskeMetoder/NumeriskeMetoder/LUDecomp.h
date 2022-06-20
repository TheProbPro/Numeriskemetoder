#ifndef LUDECOMP_H
#define LUDECOMP_H
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/utilities.h"
namespace LUDecomp {
	void LUDecomp(MatDoub& A, VecDoub& b) {
		// Evaluate x
		VecDoub x(3);
		LUdcmp LU(A);

		auto L = LU.lu;
		L[0][1] = L[0][2] = L[1][2] = 0;
		L[0][0] = L[1][1] = L[2][2] = 1;
		util::print(L, "L");

		auto U = LU.lu;
		U[1][0] = U[2][0] = U[2][1] = 0;
		util::print(U, "U");

		LU.solve(b, x);

		// print x
		util::print(A, "A");
		util::print(L * U, "L*U");
		util::print(L * U * x, "L*U*x");
		util::print(x, "x");
	}
};


#endif // LUDECOMP_H