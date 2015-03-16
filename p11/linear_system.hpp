#pragma once

#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

#include "../matrix/matrix.hpp"

/**
 * precision
 */
const double eps = 1e-5; 

/**
 * Gauss elimination
 * Does not save the matrix but modifies it in place
 * @param[in,out] mat The extended matrix of the system
 * @param[out] x Output vector containing a solution
 * @param swapVar if true, algorithm will swap variables (that is, swap columns of the matrix)
 *                when the leading element is too small
 */
template<typename T, int n>
void gauss(Math::Matrix<T, n, n+1> &mat,
		   Math::Matrix<T, n, 1> &x, bool swapVar = false) {
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
		
		for(int j = k; j < n + 1; ++j) {
			mat(k, j) /= temp;
		}

		for(int i = k + 1; i < n; ++i) {
			T temp = mat(i, k);
			for(int j = k; j < n + 1; ++j) {
				mat(i, j) -= mat(k, j) * temp;
			}
		}
	}

	/* Back traverse */
	std::array<T, n> X;
	for(int i = n - 1; i >= 0; --i) {
		T sum = 0;
		for(int j = i + 1; j < n; ++j) {
			sum += mat(i, j) * X[j];
		}
		X[i] = mat(i, n) - sum;
	}
	
	/* Restore original order */
	for(int i = 0; i < n; ++i) {
		x(newOrder[i], 0) = X[i];
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

/** Inverted matrix
 *
 */
template<typename T, int n>
void invert(Math::Matrix<T, n, n> &mat,
			Math::Matrix<T, n, n> &inv) {
	std::vector<Math::Matrix<T, n, 1> > columns;
	Math::Matrix<T, n, 1> col;

	for(int i = 0; i < n; ++i) {
		Math::Matrix<T, n, 1> identityCol;
		identityCol(i, 0) = 1;
		Math::Matrix<T, n, n+1> extended(Math::concatenateH(mat, identityCol));
		gauss(extended, col);  //LU 
		columns.push_back(col);
	}

	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < n; ++j) {
			inv(i, j) = columns[j](i, 0);
		}
	}		
}

#endif // LINEAR_SYSTEM_HPP
