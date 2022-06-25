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

		//Creating the SVD matrices
		MatDoub U = svd.u;
		MatDoub V = svd.v;
		VecDoub W = svd.w;

		//Create W_inv and prints W_inv and W
		std::cout << "Elements of W: " << std::endl;
		util::print(W);
		//Creating inverted W matrix
		MatDoub W_inv = MatDoub(W.size(), W.size());
		W_inv.assign(W.size(), W.size(), 0);
		for (int i = 0; i < W.size(); ++i) {
			W_inv[i][i] = 1. / W[i];
			if (W[i] <= 0.001) {
				W_inv[i][i] = 0;
			}
		}
		std::cout << "W inverted: " << std::endl;
		util::print(W_inv);

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

		//Create W_inv and prints W_inv and W
		std::cout << "Elements of W: " << std::endl;
		util::print(W);
		//Creating inverted W matrix
		MatDoub W_inv = MatDoub(W.size(), W.size());
		W_inv.assign(W.size(), W.size(), 0);
		for (int i = 0; i < W.size(); ++i) {
			W_inv[i][i] = 1. / W[i];
			if (W[i] <= 0.001) {
				W_inv[i][i] = 0;
			}
		}
		std::cout << "W inverted: " << std::endl;
		util::print(W_inv);

		//Solving
		ans = V * W_inv * util::T(U) * b;
		//Printing
		std::cout << "SVD solved: " << std::endl;
		util::print(ans);

		//Opgave 3 (state a vector ||y||=1 and y*(Ax) = 0 with ||x|| =1
		//----------------------Gav op!-----------------------------------------------
		/*
		std::cout << "y= " << std::endl;
		VecDoub y(10); y.assign(y.size(), 0);
		VecDoub temp(10); temp.assign(temp.size(), 0);
		for (int i = 0; i < U.ncols()-1; ++i) {
			for (int j = 0; j < U.nrows(); ++j) {
				y[j] += W_inv[i][i] * U[j][i];
			}
		}
		*/
		
	}
};

#endif // SVDDECOMP_H