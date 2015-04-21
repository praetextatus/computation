#include "../../matrix/util.hpp"

namespace Math {
	template<int n>
	double maxOverDiagonal(const Matrix<double, n, n> &A, int &idx_i, int &idx_j) {
		double max = A(0, 1);
		idx_i = 0; idx_j = 1;
		for(int i = 0; i < n-1; ++i) {
			for(int j = i+1; j < n; ++j) {
				if(std::abs(max) < std::abs(A(i, j) )) {
					max = A(i, j);
					idx_i = i;
					idx_j = j;
				}
			}
		}
		return max;
	}

	template<int n>
	void jakobi(Matrix<double, n, n> &A,
				Matrix<double, n, n> &X,
				double eps = 1e-5) {
		X = identity<double, n>();

		int i, j;
		double max = maxOverDiagonal<n>(A, i, j);

						

		while(std::abs(max) > eps) {
			auto V = identity<double, n>();
			
			double d = std::sqrt((A(i, i) - A(j, j)) * (A(i, i) - A(j, j)) + 4 * A(i, j) * A(i, j));
			double c = std::sqrt(0.5 + 0.5 * std::abs(A(i, i) - A(j, j)) / d);
			int sign = (A(i, j) * (A(i, i) - A(j, j)) > 0) ? 1 : -1;
			double s = sign * std::sqrt(0.5 - 0.5 * std::abs(A(i, i) - A(j, j)) / d);
			
			V(i, i) = c;
			V(j, j) = c;
			V(i, j) = -s;
			V(j, i) = s;

			for(int k = 0; k < n; ++k) {
				if(k == i || k == j) {
					continue;
				}

				double a_ki = A(k, i);
				double a_kj = A(k, j);

				A(k, i) = c * a_ki + s * a_kj;
				A(i, k) = A(k, i);
				A(k, j) = -s * a_ki + c * a_kj;
				A(j, k) = A(k, j);
			}
			
			double a_ii = A(i, i);
			double a_ij = A(i, j);
			double a_jj = A(j, j);
			A(i, i) = c * c * a_ii + 2 * c * s * a_ij + s * s * a_jj;
			A(j, j) = s * s * a_ii - 2 * c * s * a_ij + c * c * a_jj;
			A(i, j) = (c * c - s * s) * a_ij + c*s*(a_jj - a_ii);
			A(j, i) = A(i, j);
			X = dot(X, V);
			
			max = maxOverDiagonal<n>(A, i, j);
		}
	}	
}




