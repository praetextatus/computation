#pragma once

#include <cmath> // sqrt
#include "../matrix/matrix.hpp"
#include "../p11/linear_system.hpp"


/** A handy type.
 */
template<typename T, int n>
using Vector = Math::Matrix<T, n, 1>;


/** Euclid norm of a vector
 * @returns norm
 */
template<typename T, int n>
T euclidNorm(const Vector<T, n> &vec) {
	T sum = 0;
	for(int i = 0; i < n; ++i) {
		sum += vec(i, 0) * vec(i, 0);
	}
	return std::sqrt(sum);
}

/** Solve the system x=Hx+g iteratively.
 * @param guess The initial guess.
 * @param[out] sol Solution vector.
 */
template<typename T, int n>
int iterativeSolve(const Math::Matrix<T, n, n> &H,
				   const Vector<T, n> &g,
				   const Vector<T, n> &guess,
				   Vector<T, n> &sol,
				   double precision,
				   int maxiter = 100) {
	Vector<T, n> x = guess;
	
	for(int i = 0; i < maxiter; ++i) {
		Vector<T, n> old_x = x;
		x = Math::dot(H, x) + g;
		//std::cout << "\n" << x << "\n" << old_x;
		if(i && (euclidNorm(x - old_x) < precision)) {
			std::cout << "SHUI\n";
			sol = x;
			return i;
		}
	}
	sol = x;
	return maxiter;
}

/** Solve a system Ax=b using Seidel method.
 * Uses zero vector as the initial guess.
 * @param[out] sol Solution vector.
 */
template<typename T, int n>
int seidel(const Math::Matrix<T, n, n> &mat,
		   const Vector<T, n> &vec,
		   Vector<T, n> &sol,
		   double precision,
		   int maxiter = 100) {
	typedef Math::Matrix<T, n, n> matrix;
	matrix lower, upper, invLower;
	matrix Hseid;
	Vector<T, n> gseid;
	
	for(int row = 0; row < n; ++row) {
		for(int col = 0; col <= row; ++col) {
			lower(row, col) = mat(row, col);
		}
	}
	for(int row = 0; row < n - 1; ++row) {
		for(int col = row + 1; col < n; ++col) {
			upper(row, col) = mat(row, col);
		}
	}
	
	std::cout << "LOWER\n" << lower << "\nUPPER\n" << upper;

	invert(lower, invLower);
	Hseid = Math::dot(-invLower, upper);
	gseid = Math::dot(invLower, vec);

	Vector<T,n> guess; /* just zeros */
	
	return iterativeSolve(Hseid, gseid, guess, sol, precision, maxiter);
}
