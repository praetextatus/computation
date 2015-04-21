#include "../../matrix/linsys.hpp"
#include "../../matrix/util.hpp"

typedef Math::Matrix<double, 3, 3> Matrix3d;
typedef Math::Matrix<double, 3, 1> Vector3d;

int main() {
//	using namespace Math;
	Math::Matrix<double, 4, 4> mat 	{2, -14, 8, 1,
			3, -22, 7, 2,
			0, 2, 5, -1,
			3, -1, 3, 2};
	
	/* Invert matrix */
	std::cout << "A\n" << mat;
	Math::Matrix<double, 4, 4> inv = invert(mat);
	std::cout << "A^(-1)\n" << inv;
	
	Math::Matrix<double, 4, 1> vec{2, 0, 5, 0};
	Math::Matrix<double, 4, 1> x;
	std::cout << "\nSolving system Ax=b where b =\n" << vec;
	Math::Matrix<double, 4, 5> ext = concatenateH(mat, vec);
	gauss(ext, x);
	std::cout << "x = \n" << x;
	Math::Matrix<double, 4, 1> R = vec + (-Math::dot(mat, x));
	std::cout << "R =\n" << R;

	std::cout << "\n\nLU-Decomposition\n";
	Math::Matrix<double, 4, 4> L, U;
	luDecomposition(mat, L, U);
	std::cout << "L=\n" << L << "U\n" << U;

	mat(0,0) = mat(0,0) * 1e-8;
	std::cout << "\n\n\nSolving system Cx=b\n where C=\n" << mat;
	ext = Math::concatenateH(mat, vec);
	gauss(ext, x);	
	std::cout << "x = \n" << x;
	ext = Math::concatenateH(mat, vec);
	R = vec + (-Math::dot(mat, x));
	std::cout << "R1=\n" << R;
	gauss(ext, x, true);
	std::cout << "x = \n" << x;
	R = vec + (-Math::dot(mat, x));
	std::cout << "R2 =\n" << R;

}
	
