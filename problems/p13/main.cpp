#include <iostream>
#include "eigen.hpp"

constexpr int n = 3;

void eigen(Math::Matrix<double, n, n> mat) {
	Math::Matrix<double, n, n> eigenvec;
	Math::jakobi(mat, eigenvec, 1e-6);
	for(int i = 0; i < n; ++i) {
		std::cout << std::setprecision(6) << mat(i, i);
		std::vector<double> vec;
		for(int j = 0; j < n; ++j) vec.push_back(eigenvec(i, j));
		double norm = std::sqrt(std::accumulate(vec.begin(), vec.end(),
												0.0,
												[](double acc, double v) {
													return acc + v*v;
												}));
		std::for_each(vec.begin(), vec.end(), [norm](double &v) { v /= norm; });
		
		std::string str = std::accumulate(vec.begin()+1, vec.end(),
										  std::to_string(vec[0]),
										  [](const std::string &s, double v) {
											  return s + ", " + std::to_string(v);
										  });
		std::cout << "\t[" << str << "]\n";

	}
}
										  

int main() {
	Math::Matrix<double, 3, 3> mat {
			-0.81417, -0.01937, 0.41372,
			-0.01937, 0.54414, 0.00590,
			0.41372, 0.00590, -0.81445};
	
	Math::prettyPrint(mat, "A", 6);
	eigen(mat);
	return 0;
}
