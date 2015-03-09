#include <iostream>
#include "system2.hpp"

int main() {
	Matrix2d A{1.0, 0.99,
			0.99,  0.98};
	Vector2d b{1.99, 1.97};
	Vector2d sol = solve2(A, b);
	std::cout << "Solving system Ax=b\n" <<
		"A = \n" << A <<
		"b = \n" << b <<
		"x = \n" << sol;
	
	Vector2d b1{2.0, 2.0};
	Vector2d sol1 = solve2(A, b1);
	std::cout <<"Solving system Ax=b1\n" <<
		"A = \n" << A <<
		"b1 = \n" << b1 <<
		"x = \n" << sol1;  
	std::cout << "\nCond A = " << cond(A) << "\n";
	return 0;
}
	
	
