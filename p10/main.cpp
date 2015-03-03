#include <iostream>
#include "system2.hpp"

int main() {
	Matrix2d A{1.0, 0.99,
			0.99,  0.98};
	Vector2d b{1.99, 1.97};
	try {
		Vector2d sol = solve2(A, b);
	}
	catch(std::exception &e) {
		std::cerr << e.what();
		return 1;
	}
	std::cout << "Solving system Ax=b\n" <<
		"A = \n" << A <<
		"b = \n" << b <<
		"x = \n" << sol;
	return 0;
}
	
	
