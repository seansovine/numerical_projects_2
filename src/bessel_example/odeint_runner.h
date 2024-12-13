#ifndef ODEINT_RUNNER_H
#define ODEINT_RUNNER_H

#include <boost/numeric/odeint.hpp>

#include <vector>

// For State_T x, x[0] represents f(t) and x[1] represents f'(t).
using state_t = std::array<double, 2>;
using resultseq_t = std::vector<double>;
using timeseq_t = resultseq_t;
using results_t = std::pair<resultseq_t, timeseq_t>;

namespace detail {

// Similar to boost example, records states and times
// in the given vectors as the solver is running.

struct StateAndTimeObserver {
  std::vector<state_t> &m_states;
  std::vector<double> &m_times;

  StateAndTimeObserver(std::vector<state_t> &states, std::vector<double> &times)
      : m_states(states), m_times(times) {}

  void operator()(const state_t &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

} // namespace detail

// Concepts for OdeintRunner parameters.

template <typename T>
concept RHS = requires(T a, const state_t &x, state_t &dxdt, const double t) { a(x, dxdt, t); };

template <typename T>
concept Initializer = requires(T a) {
  { a.init() } -> std::same_as<state_t>;
};

// ------------------------------
// OdeintRunner class definition.

template <RHS Rhs, Initializer Init>
class OdeintRunner {
public:
  OdeintRunner() { x = Init{}.init(); }

  results_t run(const double inTMin, const double inTMax, const double inStep);

private:
  state_t x{};
};

// Run method definition.

template <RHS Rhs, Initializer Init>
results_t OdeintRunner<Rhs, Init>::run(const double tMin, const double tMax, const double step) {
  namespace boostode = boost::numeric::odeint;

  std::vector<state_t> x_vec;
  timeseq_t times;

  // With constant stepper and integrator, observer is called at regular intervals.
  boostode::runge_kutta4<state_t> stepper;
  boostode::integrate_const(stepper, Rhs{}, x, tMin, tMax, step, detail::StateAndTimeObserver(x_vec, times));

  // Extract function evaluations from state vec.
  resultseq_t resultsFOnly{};
  std::transform(begin(x_vec), end(x_vec), std::back_inserter(resultsFOnly),
                 [](const state_t &x) { return x[0]; });

  return {resultsFOnly, times};
}

#endif
