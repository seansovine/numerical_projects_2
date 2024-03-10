/**
 *  Compute Bessel function Jn by integrating Bessel's equation
 *  with Boost::numeric::odeint. Also compute Jn using GNU Scientific
 *  Library implementation, then plot both using matplotplusplus
 *  to compare the results.
 */

#include <utility>

#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <matplot/matplot.h>

/* Parameters. */

constexpr int BESSEL_ORDER = 4;

// Starting a little past 0 to avoid singularity.
constexpr double T_MIN = 0.01;
constexpr double T_MAX = 20.0;
constexpr double T_STEP = 0.001;

/* Definition for GSL. */

struct GslJn {
  double operator()(double x) {
    gsl_sf_result result{};
    gsl_sf_bessel_Jn_e(BESSEL_ORDER, x, &result);
    return result.val;
  }
};

/* Definitions for odeint. */

typedef std::vector<double> State_T;
typedef std::vector<double> ResultSeq_T;
typedef std::pair<ResultSeq_T, ResultSeq_T> Results_T;

struct BesselRhs {
  // RHS of equation x' = f(x).
  void operator()(const State_T &x, State_T &dxdt, const double t) {
    dxdt[0] = x[1];
    dxdt[1] = -(1 / (t * t)) * (t * x[1] + (t * t - BESSEL_ORDER * BESSEL_ORDER) * x[0]);
    // NOTE: Equation has singularity at 0.
  }
};

struct StateAndTimeObserver {
  std::vector<State_T> &m_states;
  std::vector<double> &m_times;

  StateAndTimeObserver(std::vector<State_T> &states, std::vector<double> &times)
      : m_states(states), m_times(times) {}

  void operator()(const State_T &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

class OdeintBesselRunner {
  typedef std::vector<State_T> X_Results_T_;

  // x[0] = f(t)
  // x[1] = f'(t)
  State_T x;

  // For estimating derivative.
  static constexpr double H = 0.00001;

public:
  OdeintBesselRunner() : x(2) {
    GslJn gslJn{};

    // Due to singularity in equation we're starting past x = 0, for which values are known.
    // We initialize w/ known values from GSL for verification purposes.
    x[0] = gslJn(T_MIN);
    x[1] = (gslJn(T_MIN + H) - gslJn(T_MIN - H)) / (2 * H);
  }

  Results_T run(double xMin, double xMax, double step) {
    using namespace boost::numeric::odeint;

    X_Results_T_ x_vec;
    ResultSeq_T times;

    BesselRhs bf{};

    // With constant stepper and integrator, observer is called at regular intervals.
    runge_kutta4<State_T> stepper;
    integrate_const(stepper, bf, x, xMin, xMax, step, StateAndTimeObserver(x_vec, times));

    ResultSeq_T resultsFOnly{};
    std::transform(begin(x_vec), end(x_vec), std::back_inserter(resultsFOnly),
                   [](const State_T &x) { return x[0]; });

    return std::make_pair(resultsFOnly, times);
  }
};

/* Main function. */

int main() {
  // Initialize plot.
  matplot::hold(matplot::on);
  matplot::grid(matplot::on);

  // Compute using odeint and add result to plot.

  OdeintBesselRunner runner{};
  Results_T odeintResult = runner.run(T_MIN, T_MAX, T_STEP);

  const ResultSeq_T &time = odeintResult.second;
  const ResultSeq_T &odeiVals = odeintResult.first;

  matplot::plot(time, odeiVals, "-b")->line_width(1);

  // Compute using GSL and add result to plot.

  GslJn gslJn{};
  ResultSeq_T jnVals = matplot::transform(time, gslJn);

  matplot::plot(time, jnVals, "-r")->line_width(1);

  matplot::show();
  return 0;
}
