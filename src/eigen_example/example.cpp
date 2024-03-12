
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

int main()
{
  using namespace Eigen;

  srand((unsigned int)time(0));

  DiagonalMatrix<double, 3> diag(4, 2, 1);
  Matrix3d noise = Matrix3d::Random() * 0.1;

  // Diagonal matrix perturbed by random noise
  Matrix3d m = diag.toDenseMatrix() + noise;

  // Non-perturbed eigenvector
  Vector3d v(1, 0, 0);

  // Estimate perturbed eigenvector with power method
  constexpr size_t num_iterations = 100;
  constexpr double tolerance = 1.0e-8;
  Vector3d v_last;

  for (size_t i = 0; i < num_iterations; i++) {
    v_last = v;
    v = m * v;
    v = v / v.norm();
    if ((v - v_last).norm() < tolerance) {
      break;
    }
  }

  std::cout << "Perturbed matrix:" << std::endl << m << std::endl;
  std::cout << std::endl;
  std::cout << "Perturbed eigenvector:" << std::endl << v << std::endl;
  std::cout << std::endl;
  std::cout << "Verify m * v = " << std::endl << m * v << std::endl;
  std::cout << std::endl;
  std::cout << "and " << ((m * v)(0) / v(0)) << " * v = " << std::endl
            << ((m * v)(0) / v(0)) * v << std::endl;

  return 0;
}
