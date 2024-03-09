/**
 *  Compute Bessel function Jn using GNU Scientific Library
 *  implementation, and plot using matplotplusplus.
 */

#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <matplot/matplot.h>

/* Parameters. */

constexpr int BESSEL_ORDER = 4;

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 200.0;

/* Helper definitions. */

typedef std::vector<double> State_T;
typedef std::vector<State_T> Results_T;

// Function object for RHS of x' = f(x).
class bessel {

public:
  bessel() {}

  void operator()(const State_T &x, State_T &dxdt, const double /* t */) {
    // x[0] = t, x[1] = f(t), x[2] = f'(t)
    dxdt[0] = 1;
    dxdt[1] = x[2];
    dxdt[2] = -x[1] / x[0] + (BESSEL_ORDER / (x[0] * x[0]) - 1) * x[1];
  }
};

// Observer (from odeint examples).
struct push_back_state_and_time {
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

  // state_initialization
  State_T x(3);
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;

  // Results containers.
  Results_T x_vec;
  std::vector<double> times;

  bessel bf{};

  // define_adapt_stepper
  typedef runge_kutta_cash_karp54<State_T> error_stepper_type;

  // integrate_adapt
  typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;
  integrate_adaptive(controlled_stepper, bf, x, xMin, xMax, step, push_back_state_and_time(x_vec, times));

  // TOOD: Transform it to just return f values.
  return x_vec;
}

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
