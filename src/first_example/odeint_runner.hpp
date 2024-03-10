#ifndef ODEINT_RUNNER_H
#define ODEINT_RUNNER_H

#include <boost/numeric/odeint.hpp>

#include <vector>

// For State_T x, x[0] represents f(t) and x[1] represents f'(t).
typedef std::vector<double> State_T;
typedef std::vector<double> ResultSeq_T;
typedef std::pair<ResultSeq_T, ResultSeq_T> Results_T;

template <class Rhs, class Initializer>
class OdeintRunner {

public:
  OdeintRunner(double inTMin, double inTMax, double inStep) : tMin{inTMin}, tMax{inTMax}, step{inStep}, x(2) {
    // Since we're starting past x = 0, where values are known, we
    // initialize w/ known values from GSL for verification purposes.

    Initializer initFunc{};

    x[0] = initFunc(tMin);
    x[1] = (initFunc(tMin + H) - initFunc(tMin - H)) / (2 * H);
  }

  Results_T run();

private:
  typedef std::vector<State_T> X_Results_T_;

  struct StateAndTimeObserver;

  const double tMin;
  const double tMax;
  const double step;

  State_T x;

  // For estimating derivative.
  static constexpr double H = 0.00001;
};

template <class Rhs, class Initializer>
struct OdeintRunner<Rhs, Initializer>::StateAndTimeObserver {
  std::vector<State_T> &m_states;
  std::vector<double> &m_times;

  StateAndTimeObserver(std::vector<State_T> &states, std::vector<double> &times)
      : m_states(states), m_times(times) {}

  void operator()(const State_T &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

template <class Rhs, class Initializer>
Results_T OdeintRunner<Rhs, Initializer>::run() {
  using namespace boost::numeric::odeint;

  X_Results_T_ x_vec;
  ResultSeq_T times;

  // With constant stepper and integrator, observer is called at regular intervals.
  runge_kutta4<State_T> stepper;
  integrate_const(stepper, Rhs{}, x, tMin, tMax, step, StateAndTimeObserver(x_vec, times));

  ResultSeq_T resultsFOnly{};
  std::transform(begin(x_vec), end(x_vec), std::back_inserter(resultsFOnly),
                 [](const State_T &x) { return x[0]; });

  return std::make_pair(resultsFOnly, times);
}

#endif
