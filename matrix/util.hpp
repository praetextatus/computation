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
};
