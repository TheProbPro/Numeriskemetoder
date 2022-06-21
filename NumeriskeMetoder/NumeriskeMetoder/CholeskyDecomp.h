#ifndef CHOLESKYDECOMP_H
#define CHOLESKYDECOMP_H

//Includes
#include "../../NR_C301/code/nr3.h"
//#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/cholesky.h"
#include "../../NR_C301/code/utilities.h"

namespace CholeskyDecomp
{
	void CholeskyDecomp(MatDoub &C, VecDoub &c, int size){
		// Createing vector to contain answer
		VecDoub chCoef(c.size());
		// Solving cholesky
		Cholesky ch(C);
		ch.solve(c, chCoef);

		// Printing answer
		std::cout << "Cholesky Decomposition answers: " << std::endl;
		util::print(chCoef);
	}
};

#endif // CHOLESKYDECOMP_H