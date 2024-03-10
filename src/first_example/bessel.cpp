/**
 *  Compute Bessel function Jn by integrating Bessel's equation
 *  with Boost::numeric::odeint. Also compute Jn using GNU Scientific
 *  Library implementation, then plot both using matplotplusplus
 *  to compare the results.
 */

#include <utility>

#include <gsl/gsl_sf_bessel.h>
#include <matplot/matplot.h>

#include "odeint_runner.hpp"

// Parameters.

constexpr int BESSEL_ORDER = 4;

// Start past 0 to avoid singularity.
constexpr double T_MIN = 0.01;
constexpr double T_MAX = 20.0;
constexpr double T_STEP = 0.001;

// GSL helper.

struct GslBesseln {
  double operator()(double x) {
    gsl_sf_result result{};
    gsl_sf_bessel_Jn_e(BESSEL_ORDER, x, &result);
    return result.val;
  }
};

// Odeint RHS.

struct BesselRhs {
  // RHS of equation x' = f(x).
  void operator()(const State_T &x, State_T &dxdt, const double t) {
    dxdt[0] = x[1];
    dxdt[1] = -(1 / (t * t)) * (t * x[1] + (t * t - BESSEL_ORDER * BESSEL_ORDER) * x[0]);
    // NOTE: Equation has singularity at 0.
  }
};

using OdeintBesselRunner = OdeintRunner<BesselRhs, GslBesseln>;

// Matplot helper.

struct MatplotStateManager {
  MatplotStateManager() {
    // Initialize plot.
    matplot::hold(matplot::on);
    matplot::grid(matplot::on);
  }

  ~MatplotStateManager() { matplot::show(); }
};

// Main.

int main() {
  Results_T odeintResult = OdeintBesselRunner(T_MIN, T_MAX, T_STEP).run();
  const ResultSeq_T &time = odeintResult.second;
  const ResultSeq_T &odeintVals = odeintResult.first;

  ResultSeq_T gslJnVals = matplot::transform(time, GslBesseln{});

  {
    // Configures plot and shows on scope exit.
    MatplotStateManager mgr{};

    matplot::plot(time, odeintVals, "-b")->line_width(1);
    matplot::plot(time, gslJnVals, "-r")->line_width(1);
  }

  return 0;
}
