#include <vector>

#include <boost/numeric/odeint.hpp>

// rhs_function
/* The type of container used to hold the state vector */
typedef std::vector<double> state_type;

const int n = 4;

/* The rhs of x' = f(x) defined as a class */
class bessel {

public:
  bessel() {}

  void operator()(const state_type &x, state_type &dxdt, const double /* t */) {
    // x[0] = t
    // x[1] = f(t)
    // x[2] = f'(t)

    dxdt[0] = 1;
    dxdt[1] = x[2];
    dxdt[2] = -x[1] / x[0] + (n / (x[0] * x[0]) - 1) * x[1];
  }
};

int main() {
  using namespace boost::numeric::odeint;

  bessel bf{};

  // state_initialization
  state_type x(2);
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;

  // define_adapt_stepper
  typedef runge_kutta_cash_karp54<state_type> error_stepper_type;

  // integrate_adapt
  typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;
  integrate_adaptive(controlled_stepper, bf, x, 0.0, 10.0, 0.01);
}
