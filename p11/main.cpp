#include <iostream>
#include "linear_system.hpp"

typedef Math::Matrix<double, 3, 3> Matrix3d;
typedef Math::Matrix<double, 3, 1> Vector3d;

int main() {
	Matrix3d mat 	{2, -14, 8,
			3, -22, 7,
			0, 2, 5};
	
	/* Invert matrix */
	std::cout << "A\n" << mat;
	Matrix3d inv;
	invert(mat, inv);
	std::cout << "A^(-1)\n" << inv;
	
	Vector3d vec{2, 0, 5};
	Vector3d x;
	std::cout << "\nSolving system Ax=b where b =\n" << vec;
	Math::Matrix<double, 3, 4> ext = Math::concatenateH(mat, vec);
	gauss(ext, x);
	std::cout << "x = \n" << x;
	Vector3d R = vec + (-Math::dot(mat, x));
	std::cout << "R =\n" << R;

	std::cout << "\n\nLU-Decomposition\n";
	Math::Matrix<double, 3, 3> L, U;
	luDecomposition(mat, L, U);
	std::cout << "L=\n" << L << "U\n" << U;

	mat(0,0) = mat(0,0) * 1e-8;
	std::cout << "\n\n\nSolving system Cx=b\n where C=\n" << mat;
	ext = Math::concatenateH(mat, vec);
	gauss(ext, x, true);
	std::cout << "x = \n" << x;
	R = vec + (-Math::dot(mat, x));
	std::cout << "R =\n" << R;

}
	
