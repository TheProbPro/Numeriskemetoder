#ifndef SVDDECOMP_H
#define SVDDECOMP_H

// Includes
#include "../../NR_C301/code/nr3.h"
#include "../../NR_C301/code/utilities.h"
#include "../../NR_C301/code/svd.h"

namespace SVDDecomp
{
	void SVDDecomp(MatDoub& A, VecDoub& b, int nCols) {
		//Create SVD object and answer vector
		SVD svd(A);
		VecDoub ans(nCols);
		//Solving and printing
		svd.solve(b, ans, svd.eps);

		std::cout << "SVD solved: " << std::endl;
		util::print(ans);
	}

	void SVDSolveAlg(MatDoub& A, VecDoub& b, int nCols) {
		//Create SVD object and answer vector
		SVD svd(A);
		VecDoub ans(nCols);

		//Writing the solver
		MatDoub U = svd.u;
		MatDoub V = svd.v;
		VecDoub W = svd.w;

		//Creating inverted W matrix
		MatDoub Winv = MatDoub(W.size(), W.size());
		for (int i = 0; i < W.size(); i++) {
			Winv[i][i] = 1.0 / W[i];
		}

		//Solving
		ans = V * Winv * util::T(U) * b;
		//Printing
		std::cout << "SVD solved: " << std::endl;
		util::print(ans);
	}
};

#endif // SVDDECOMP_H