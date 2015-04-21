#include <iostream>
#include "../../matrix/linsys.hpp"
#include "../../matrix/util.hpp"

constexpr int n = 3;
typedef Math::Matrix<double, n, n> Matrix;
typedef Math::Matrix<double, n, 1> Vector;

int main() {
	Matrix mat 	{4, -1, -1,
			-2, 6, 1,
			-1, 1, 7};
	
	Math::prettyPrint(mat, "A");

    Vector vec{3, 9, -6};
	Vector x;
	std::cout << "\nSolving system Ax=b with Gauss where\n";
	Math::prettyPrint(vec, "b");;
	Math::Matrix<double, 3, 4> ext = Math::concatenateH(mat, vec);
	gauss(ext, x);
	Math::prettyPrint(x, "x");
	Vector R = vec -Math::dot(mat, x);
	Math::prettyPrint(R, "R");

	std::cout << "\nSolving iteratively\n";
	x = {0, 0, 0};
	Matrix H;
	Vector g;
	rewriteSystem(mat, vec, H, g);
	int iter = iterativeSolve(H, g, {0, 0, 0}, x);
	Math::prettyPrint(x, "x");
	R = vec - Math::dot(mat, x);
	Math::prettyPrint(R, "R");
	std::cout << "iterations = " << iter << "\n";


	x = {0, 0, 0};
	std::cout << "\nSolving the same system with Seidel\n";
	iter = seidel(mat, vec, x, 1e-5);
	Math::prettyPrint(x, "x");
	R = vec -Math::dot(mat, x);
	Math::prettyPrint(R, "R");
	std::cout << "iterations = " << iter << "\n";

}
	
