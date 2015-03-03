#include <boost/test/unit_test.hpp>
#include "../linear_system.hpp"

BOOST_AUTO_TEST_SUITE( test_suite_linsys );

BOOST_AUTO_TEST_CASE( test_gauss ) {
	/* 2x + 3y = 6
	   4x + 9y = 15 */
	Math::Matrix<double, 2, 3> 
		mat {2, 3,  6,
			4, 9, 15}; 
	Math::Matrix<double, 2, 1> expected_x{1.5, 1};
	Math::Matrix<double, 2, 1> x;
	
	std::cout << "SOLVING\n";
	std::cout << mat;
	std::cout << "EXPECTED\n" << expected_x;
	gauss(mat, x);
	std::cout << "GOT\n" << x;
	std::cout << "GOT\n" << mat;
	BOOST_CHECK(x == expected_x);
}

BOOST_AUTO_TEST_CASE( test_LU ) {
	Math::Matrix<double, 3, 3> 
		mat {2, -14, 8,
			3, -22, 7,
			0, 2, 5};
	Math::Matrix<double, 3, 3>
		L_exp {2, 0, 0,
			3, -1, 0,
			0, 2, -5};
	Math::Matrix<double, 3, 3>
		U_exp {1, -7, 4, 
			0, 1, 5,
			0, 0, 1};
	Math::Matrix<double, 3, 3> L, U;
	std::cout << "LU-decomposition\n";
	std::cout << mat;
	std::cout << "EXPECTED\n";
	std::cout << "L\n" << L_exp;
	std::cout << "U\n" << U_exp;
	luDecomposition(mat, L, U);
	std::cout << "GOT\n";
	std::cout << "L\n" << L;
	std::cout << "U\n" << U;
	BOOST_CHECK(L == L_exp);
	BOOST_CHECK(U == U_exp);
}
	

BOOST_AUTO_TEST_SUITE_END();
