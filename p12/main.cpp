#include <iostream>
#include "seidel.hpp"

constexpr int n = 3;
typedef Math::Matrix<double, n, n> Matrix3;
typedef Math::Matrix<double, n, 1> Vector3;

int main() {
	Matrix3 mat 	{4, -1, -1,
			-2, 6, 1,
			-1, 1, 7};
	
	std::cout << mat; 

    Vector3 vec{3, 9, -6};
	Vector3 x;
	std::cout << "\nSolving system Ax=b with Gauss where b =\n" << vec;
	Math::Matrix<double, 3, 4> ext = Math::concatenateH(mat, vec);
	gauss(ext, x);
	std::cout << "x = \n" << x;
	Vector3 R = vec -Math::dot(mat, x);
	std::cout << "R =\n" << R;

	x = {0, 0, 0};
	std::cout << "\nSolving the same system with Seidel\n";
	int iter = seidel(mat, vec, x, 1e-5);
	std::cout << "x = \n" << x;
	std::cout << "iterations = " << iter << "\n";
	R = vec -Math::dot(mat, x);
	std::cout << "R =\n" << R;
}
	
