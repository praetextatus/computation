#pragma once

#include "linsys.hpp"

namespace Math {
	/**
	 * @returns Identity matrix.
	 */
	template<typename T, int n>
	Math::Matrix<T, n, n> identity() {
		Math::Matrix<T, n, n> E;
		for(int i = 0; i < n; ++i) {
			E(i, i) = 1;
		}
		return E;
	}

	
	/** Inverted matrix
	 *
	 */
	template<typename T, int n>
	Math::Matrix<T, n, n> invert(const Math::Matrix<T, n, n> &mat) {
		Math::Matrix<T, n, n> inv;
		Math::Matrix<T, n, n+n> cat = concatenateH(mat, identity<T, n>());
		gauss(cat, inv);
		return inv;
	}

	/**
	 * Concatenate two matrices horizontally.
	 * @param mat1 Left-hand nxr matrix.
	 * @param mat2 Right-hand nxs matrix.
	 * @returns A nx(r+s) matrix of which first r columns are from mat1 and last s columns are from mat2.
	 */
	template<typename T, int n, int r, int s>
	Matrix<T, n, r+s> concatenateH(const Matrix<T, n, r> &mat1, const Matrix<T, n, s> &mat2) {
		Matrix<T, n, r+s> cat;
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < r + s; ++j) {
				cat(i, j) = (j < r) ? mat1(i, j) : mat2(i, j-r);
			}
		}
		return cat;
	}

	/**
	 * Pretty print
	 * Prints matrix in a pretty way. New line is inserted after.
	 */
	template<typename T, int n, int m>
	std::ostream& operator<<(std::ostream &os, const Matrix<T, n, m> &mat) {
		os << std::setprecision(3);
		for(int i = 0; i < n; ++i) {
			os << "| ";
			for(int j = 0; j < m; ++j) {
				if(std::abs(mat(i, j)) < 0.1 && std::abs(mat(i, j)) > 1e6) {
					os << std::setiosflags(std::ios::scientific);
				}
				os << std::setw(10) 
				   << mat(i, j) << " ";
			}
			os << "|\n";
		}
		return os;
	}

	/**
	 * Pretty print
	 * A better pretty print. Writes directly to std::cout. New line is inserted after.
	 */
	template<typename T, int n, int m>
    void prettyPrint(const Matrix<T, n, m> &mat, std::string label = "", int precision = 3) {
		std::cout << std::setprecision(precision);
		int labelLen = label.length();
		for(int i = 0; i < n; ++i) {
			if(labelLen) {
				std::cout << ((i == (n-1)/2)
							  ? (label + " = ")
							  : std::string(labelLen + 3, ' '));
			}
			std::cout << "| ";
			for(int j = 0; j < m; ++j) {
				std::cout << std::setw(10);
				std::cout << mat(i, j) << " ";
			}
			std::cout << "|\n";
		}
		std::cout << "\n";
	}
};
