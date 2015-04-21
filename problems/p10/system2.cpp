#include "system2.hpp"



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

double norm2(const Matrix2d &mat) {
	double norm = 0;
	for(int i = 0; i < 2; ++i) {
		double sum = 0;
		for(int j = 0; j < 2; ++j) {
			sum += std::abs(mat(i,j));
		}
		if(sum > norm) {
			norm = sum;
		}
	}
	return norm;
}

double cond(const Matrix2d &mat) {
	Matrix2d inv = inverse2(mat);
	return norm2(mat) * norm2(inv);
}
