#pragma once

#include <cmath> // sqrt
#include "../matrix/matrix.hpp"
#include "../p11/linear_system.hpp"


/** A handy type.
 */
//template<typename T, int n>
//using Vector = Math::Matrix<T, n, 1>;


/** Euclid norm of a vector
 * @returns norm
 */
template<typename T, int n>
T euclidNorm(const Math::Matrix<T, n, 1> &vec) {
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
				   const Math::Matrix<T, n, 1> &g,
				   const Math::Matrix<T, n, 1> &guess,
				   Math::Matrix<T, n, 1> &sol,
				   double precision,
				   int maxiter = 100) {
	Math::Matrix<T, n, 1> x = guess;
	
	for(int i = 0; i < maxiter; ++i) {
		Math::Matrix<T, n, 1> old_x = x;
		x = Math::dot(H, x) + g;
		//std::cout << "\n" << x << "\n" << old_x;
		if(i && (euclidNorm(x - old_x) < precision)) {
			sol = x;
			return i;
		}
	}
	sol = x;
	return maxiter;
}

template<typename T, int n>
void rewriteSystem(const Math::Matrix<T, n, n> &mat,
				   const Math::Matrix<T, n, 1> &b,
				   Math::Matrix<T, n, n> &H,
				   Math::Matrix<T, n, 1> &g) {
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < n; ++j) {
			H(i, j) = (i == j) ? 0 : -mat(i, j) / mat(i, i);
			g(i, 0) = b(i, 0) / mat(i, i);
		}
	}
}

/** Solve a system \p mat*sol=vec using Seidel method.
 * Uses zero vector as the initial guess.
 * @param[out] sol Solution vector.
 */
template<typename T, int n>
int seidel(const Math::Matrix<T, n, n> &mat,
		   const Math::Matrix<T, n, 1> &vec,
		   Math::Matrix<T, n, 1> &sol,
		   double precision,
		   int maxiter = 100) {
	typedef Math::Matrix<T, n, n> matrix;
	typedef Math::Matrix<T, n, 1> vector;

	matrix Hseid;
	vector gseid;

	
	rewriteSystem(mat, vec, Hseid, gseid);

	std::cout << Hseid << "\n" << gseid;
	
	vector x;
	vector old_x;
	for(int iter = 0; iter < maxiter; ++iter) {
		old_x = x;
		for(int i = 0; i < n; ++i) {
			T sum = 0;
			for(int j = 0; j < n; ++j) {
				sum += (j < i) ? Hseid(i, j) * x(j, 0) : Hseid(i, j) * old_x(j, 0);
			}
			sum += gseid(i, 0);
			x(i, 0) = sum;
		}
		if(iter && (euclidNorm(x - old_x) < eps)) {
			sol = x;
			return iter;
		}
	}
	
	sol = x;
	return maxiter;
}
