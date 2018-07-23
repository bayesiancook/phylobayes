
#ifndef LINALG_H
#define LINALG_H

class LinAlg {

	public:

	// diagonalize a reversible rate matrix
	// first transforms reversible matrix into a symmetric matrix
	// then calls DiagonalizeSymmetricMatrix
	static int DiagonalizeRateMatrix(double** u, double* pi, int dim, double* eigenval, double** eigenvect, double** inveigenvect, int nmax=1000, double epsilon = 1e-10);

	// diagonalize a symmetric matrix
	// first applying Householder transformation (tri-diagonal)
	// then using QR reduction

	// nmax : max number of iterations
	// epsilon : loops until non-diagonal elements are < epsilon in absolute value
	// returns number of iterations
	// eigenvect matrix is orthonormal (its transpose is its inverse)
	static int DiagonalizeSymmetricMatrix(double** u, int dim, int nmax, double epsilon, double* eigenval, double** eigenvect);

	// computes inverse of matrix given as an input (a)
	// by Gauss elimination
	// store inverse in invu
	// does not corrupt matrix a
	// if invu not specified, just returns the logdet
	static double Gauss(double** a, int dim, double** invu = 0);

	private:

	static void QR(double** u, int dim, double** ql, double** r);
	static void HouseHolder(double** u, int dim, double** a, double** ql);

};

#endif
