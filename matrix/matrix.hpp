#pragma once

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <initializer_list>
#include <array>
#include <assert.h>
#include <iomanip>

namespace Math 
{
	/**
	 * A matrix class
	 * @tparam T type of matrix elements
	 * @tparam n number of rows
	 * @tparam m number of columns
	 */
	template<typename T, int n, int m>
	class Matrix 
	{
	protected:
		std::array<T,n*m> matrix;
	public:
		Matrix();
		Matrix(std::initializer_list<T> l);
		Matrix(const Matrix<T, n, m> &mat);
		Matrix(Matrix<T, n, m> &&mat);
		T& operator()(int row, int col);
		T operator()(int row, int col) const;
		Matrix<T, n, m>& operator=(Matrix<T, n, m> mat);
		friend void swap(Matrix& m1, Matrix& m2) {
			std::swap(m1.matrix, m2.matrix);
		};

		/* arithmetics */
		Matrix<T, n, m>& operator+=(const Matrix<T, n, m>& rhs);
		Matrix<T, n, m>& operator*=(T scalar);
		bool operator==(const Matrix<T, n, m> &rhs);
	};

	template<typename T, int n, int m>
	Matrix<T, n, m>::Matrix() 
	{
		assert(n > 0 && m > 0);
		matrix.fill(0);
	}

	template<typename T, int n, int m>
	Matrix<T, n, m>::Matrix(std::initializer_list<T> l) 
	{
		assert(l.size() <= n*m);
		std::copy(l.begin(), l.end(), matrix.begin());
	}

	template<typename T, int n, int m>
	Matrix<T, n, m>::Matrix(const Matrix<T, n, m> &mat) {
		std::copy(mat.matrix.begin(), mat.matrix.end(), matrix.begin());
	}

	template<typename T, int n, int m>
	Matrix<T, n, m>::Matrix(Matrix<T, n, m> &&mat) 
	{
		swap(*this, mat);
	}
		

	template<typename T, int n, int m>
	T& Matrix<T, n, m>::operator()(int row, int col) 
	{
		assert(row >= 0 && row < n);
		assert(col >= 0 && col < m);
		return matrix[row*m+col];
	}

	template<typename T, int n, int m>
	T Matrix<T, n, m>::operator()(int row, int col) const {
		assert(row >= 0 && row < n);
		assert(col >= 0 && col < m);
		return matrix[row*m+col];
	}

	template<typename T, int n, int m>
	Matrix<T, n, m>& Matrix<T, n, m>::operator=(Matrix<T, n, m> mat) {
		swap(*this, mat);
		return *this;
	}

	/* arithmetics */
	template<typename T, int n, int m>
	Matrix<T, n, m>& Matrix<T, n, m>::operator+=(const Matrix<T, n, m>& rhs) {
		for(int i = 0; i < n*m; ++i) {
			matrix[i] += rhs.matrix[i];
		}
		return *this;
	}
	
	template<typename T, int n, int m>
	Matrix<T, n, m>& Matrix<T, n, m>::operator*=(T scalar) {
		for(auto it = matrix.begin(); it != matrix.end(); ++it) {
			(*it) *= scalar;
		}
		return *this;
	}
	
	template<typename T, int n, int m>
	Matrix<T, n, m> operator-(const Matrix<T, n, m> &rhs) {
		Matrix<T, n, m> result;
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < m; ++j) {
				result(i, j) = -rhs(i, j);
			}
		}
		return result;
	}

	template<typename T, int n, int m>
	bool Matrix<T,n,m>::operator==(const Matrix<T, n, m> &rhs) {
		bool result = true;
		for(int i = 0; i < n*m; ++i) {
			if(matrix[i] != rhs.matrix[i]) {
				result = false;
				break;
			}
		}
		return result;
	}

	template<typename T, int n, int m>
	Matrix<T, n, m> operator+(Matrix<T, n, m> lhs, const Matrix<T, n, m> &rhs) {
		lhs += rhs;
		return lhs;
	}

	/**
	 * A scalar multiplication.
	 * @see dot()
	 */
	template<typename T, int n, int m>
	Matrix<T, n, m> operator*(Matrix<T, n, m> lhs, T scalar) {
		lhs *= scalar;
		return lhs;
	}

	template<typename T, int n, int m>
	Matrix<T, n, m> operator*(T scalar, const Matrix<T, n, m> &rhs) {
		return rhs * scalar;
	}

	/**
	 * Concatenate two matrices horizontally.
	 * That is, given a matrix (A) and a matrix (B) a matrix (A B) will be returned.
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
	 * Dot product.
	 * This is a naive algorithm. It will be replaced by some other algorithm in the future 
	 */
	template<typename T, int n, int m, int r>
	Matrix<T, n, m> dot(const Matrix<T, n, r> &lhs, const Matrix<T, r, m> &rhs) {
		Matrix<T, n, m> product;
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < m; ++j) {
				for (int k = 0; k < r; ++k) {
					product(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		return product;
	}

	/**
	 * Output operator.
	 * Prints matrix in a pretty way. New line is inserted after.
	 */
	template<typename T, int n, int m>
	std::ostream& operator<<(std::ostream &os, const Matrix<T, n, m> &mat) {
		for(int i = 0; i < n; ++i) {
			os << "| ";
			for(int j = 0; j < m; ++j) {
				os << std::setprecision(5) 
				   << std::setw(7)
				   << mat(i, j) << " ";
			}
			os << "|\n";
		}
		return os;
	}
};

#endif // MATRIX_HPP
