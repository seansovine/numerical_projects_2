#include <gsl/gsl_sf_bessel.h>
#include <matplot/matplot.h>

constexpr int BESSEL_ORDER = 2;

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 500.0;

int main() {
  auto Jn = [](double x) {
    gsl_sf_result result{};
    gsl_sf_bessel_Jn_e(BESSEL_ORDER, x, &result);
    return result.val;
  };

  matplot::hold(matplot::on);

  matplot::fplot(Jn, std::array<double, 2>{X_MIN, X_MAX}, "-r")->line_width(1);
  matplot::grid(matplot::on);

  matplot::show();
  return 0;
}
