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

// Start a little past 0 to avoid singularity.
constexpr double T_MIN = 2;
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

class Bessel {
  // Function object for RHS of x' = f(x).

public:
  Bessel() {}

  void operator()(const State_T &x, State_T &dxdt, const double t) {
    dxdt[0] = x[1];
    dxdt[1] = -(1 / (t * t)) * (t * x[1] + (t * t - BESSEL_ORDER * BESSEL_ORDER) * x[0]);
    // NOTE: Singularity at 0.
  }
};

struct push_back_state_and_time {
  // Observer (from odeint examples).

  std::vector<State_T> &m_states;
  std::vector<double> &m_times;

  push_back_state_and_time(std::vector<State_T> &states, std::vector<double> &times)
      : m_states(states), m_times(times) {}

  void operator()(const State_T &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

Results_T runOdeint(double xMin, double xMax, double step) {
  using namespace boost::numeric::odeint;
  typedef std::vector<State_T> X_Results_T_;

  GslJn gslJn{};

  constexpr double H = 0.00001;

  // state_initialization
  State_T x(2);
  x[0] = gslJn(T_MIN);
  x[1] = (gslJn(T_MIN + H) - gslJn(T_MIN - H)) / (2 * H);
  // Initialize w/ known value for testing.

  // Results containers.
  X_Results_T_ x_vec;
  ResultSeq_T times;

  Bessel bf{};

  // Use a simple constant stepper and integrator.
  // This is because we want the observer to be called
  // at regular intervals. (See odeint docs.)
  runge_kutta4<State_T> stepper;
  integrate_const(stepper, bf, x, xMin, xMax, step, push_back_state_and_time(x_vec, times));

  ResultSeq_T resultsFOnly(size(times));
  for (size_t i = 0; i < size(times); i++) {
    resultsFOnly[i] = x_vec[i][0];
  }

  return std::make_pair(resultsFOnly, times);
}

/* Main function. */

int main() {
  matplot::hold(matplot::on);
  matplot::grid(matplot::on);

  // Compute and plot using odeint.

  Results_T odeintResult = runOdeint(T_MIN, T_MAX, T_STEP);

  const ResultSeq_T &time = odeintResult.second;
  const ResultSeq_T &odeiVals = odeintResult.first;

  matplot::plot(time, odeiVals, "-b")->line_width(1);

  // Compute and plot using GSL.

  GslJn gslJn{};

  ResultSeq_T jnVals = matplot::transform(time, gslJn);

  matplot::plot(time, jnVals, "-r")->line_width(1);

  // Show the plot.

  matplot::show();
  return 0;
}
