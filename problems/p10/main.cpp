#include <iostream>
#include "system2.hpp"

double norm1(const Vector2d &v) {
	return std::abs(v(0, 0)) + std::abs(v(1,0));
}

int main() {
	Matrix2d A{1.0, 0.99,
			0.99,  0.98};
	Vector2d b{1.99, 1.97};
	Vector2d sol = solve2(A, b);
	std::cout << "Solving system Ax=b\n" <<
		"A = \n" << A <<
		"b = \n" << b <<
		"x = \n" << sol;

	std::cout << "A^(-1) =\n" << inverse2(A) << "\n";

	std::cout << "Norm A = " << norm2(A) << "\n";
	std::cout << "Norm A^(-1) = " << norm2(inverse2(A)) << "\n";
	
	Vector2d b1{2.0, 2.0};
	Vector2d sol1 = solve2(A, b1);
	std::cout <<"Solving system Ax=b1\n" <<
		"A = \n" << A <<
		"b1 = \n" << b1 <<
		"x = \n" << sol1;  
	double condA = cond(A);
	std::cout << "\nCond A = " << condA << "\n";

	Vector2d dsol = sol1 + (-sol);
	Vector2d db = b1 + (-b);
	double err = norm1(dsol) / norm1(sol);
	double estim = condA * (norm1(db)/ norm1(b));
	std::cout << err << " <= " << estim << "\n";

	return 0;
}
	
	
