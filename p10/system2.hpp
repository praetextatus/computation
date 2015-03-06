#pragma once

#ifndef SYSTEM2_HPP
#define SYSTEM2_HPP

#include <cmath>
#include <stdexcept>
#include "../matrix/matrix.hpp"

typedef Math::Matrix<double, 2, 2> Matrix2d;
typedef Math::Matrix<double, 2, 1> Vector2d;

/**
 * Determinant of a 2x2 matrix 
 */
inline double det2(const Matrix2d &mat) {
	return mat(0,0) * mat(1,1) - mat(0,1) * mat(0,1);
}

/**
 * Inverse of a 2x2 matrix
 */
Matrix2d inverse2(const Matrix2d &mat);

/**
 * Solve 2x2 system of linear equations
 */
Vector2d solve2(const Matrix2d &coeffs, const Vector2d &b);

/**
 * Norm of a matrix
 */
double norm2(const Matrix2d &mat);


/**
 * Conditional number
 */
double cond(const Matrix2d &mat);
	
	
	
#endif // SYSTEM2_HPP	
