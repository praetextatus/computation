#include "system2.hpp"

double det2(const Matrix2d &mat) {
	return mat(0,0) * mat(1,1) - mat(0,1) * mat(0,1);
}

Matrix2d inverse2(const Matrix2d &mat) {
	double det = det2(mat);
	if(det == 0) {
		throw std::domain_error("Not invertible matrix");
	}
	Matrix2d inv{ mat(1,1), -mat(1,0), 
			-mat(0,1), mat(0,0)};
	inv *= 1.0/det;
	return inv;
}

Vector2d solve2(const Matrix2d &coeffs, const Vector2d &b){
	Matrix2d inv = inverse2(coeffs);
	Vector2d solution = Math::dot(inv, b);
	return solution;
}
