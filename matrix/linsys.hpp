#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

#include "matrix.hpp"

namespace Math {
	/**
	 * Gauss elimination.
	 * Does not save the matrix but modifies it in place.
	 * It can be used to solve `s` linear systems A*x_i = b_i (i = 0..s-1) simultaneously.
	 * @param[in,out] mat The extended matrix of the systems. First n rows are A, last s rows are b_i.
	 * @param[out] x Output matrix containing solutions x_i.
	 * @param swapVar if true, algorithm will swap variables (that is, swap columns of the matrix)
	 *                when the leading element is too small.
	 * @param eps precision
	 */
	template<typename T, int n, int s>
	void gauss(Math::Matrix<T, n, n+s> &mat,
			   Math::Matrix<T, n, s> &x, bool swapVar = false, double eps = 1e-5) {
		std::array<int, n> newOrder;
		for(int i = 0; i < n; ++i) {
			newOrder[i] = i;
		}
	
		/* Direct traverse */
		for(int k = 0; k < n; ++k) {
			T temp = mat(k, k);
			if(std::abs(temp) < eps) {
				std::cout << boost::format("Small leading element %1$e\n") % temp;
			}
			if(swapVar) {
				T curMax = temp;
				int maxIdx = k;
				for(int i = k + 1; i < n; ++i) {
					if(std::abs(mat(k, i)) > std::abs(curMax)) {
						curMax = mat(k, i);
						maxIdx = i;
					}
				}
				if(maxIdx != k) {
					for(int z = 0; z < n; ++z) {
						std::swap(mat(z,k), mat(z, maxIdx));
					}
					temp = mat(k, k);
					std::swap(newOrder[k], newOrder[maxIdx]);
					std::cout << "Swapped columns\n";
				}
			}			
		
			for(int j = k; j < n + s; ++j) {
				mat(k, j) /= temp;
			}

			for(int i = k + 1; i < n; ++i) {
				T temp = mat(i, k);
				for(int j = k; j < n + s; ++j) {
					mat(i, j) -= mat(k, j) * temp;
				}
			}
		}

		/* Back traverse */
		Math::Matrix<T, n, s> X;
		for(int k = 0; k < s; ++k) {
			for(int i = n - 1; i >= 0; --i) {
				T sum = 0;
				for(int j = i + 1; j < n; ++j) {
					sum += mat(i, j) * X(j,k);
				}
				X(i,k) = mat(i, n+k) - sum;
			}
		}
	
		/* Restore original order */
		for(int k = 0; k < s; ++k) {
			for(int i = 0; i < n; ++i) {
				x(newOrder[i], k) = X(i, k);
			}
		}
	}

	/** LU-decomposition
	 * @param[in] mat Matrix to be decomposed
	 * @param[out] L lower-triangular matrix; no clearing is done (TODO?)
	 * @param[out] U upper-unitriangular matrix; no clearing is done (TODO?)
	 */
	template<typename T, int n>
	void luDecomposition(const Math::Matrix<T, n, n> &mat,
						 Math::Matrix<T, n, n> &L,
						 Math::Matrix<T, n, n> &U) {
		for(int i = 0; i < n; ++i) {
			for(int j = i; j < n; ++j) {
				T sum = 0;
				for(int k = 0; k < j; ++k) {
					sum += L(j, k) * U(k, i);
				}
				L(j, i) = mat(j, i) - sum;
				sum = 0;
				for(int k = 0; k < i; ++k) {
					sum += L(i, k) * U(k, j);
				}
				U(i, j) = (mat(i, j) - sum) / L(i, i);
			}
		}
	}

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
					   double precision = 1e-5,
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

	/** Rewrite system Ax=b to form x=Hx+g */
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
			   double eps = 1e-5,
			   int maxiter = 100) {
		typedef Math::Matrix<T, n, n> matrix;
		typedef Math::Matrix<T, n, 1> vector;

		matrix Hseid;
		vector gseid;

		rewriteSystem(mat, vec, Hseid, gseid);

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

};
