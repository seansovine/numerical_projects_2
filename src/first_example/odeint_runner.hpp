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
  OdeintRunner() : x(2) { x = Initializer{}.init(); }

  Results_T run(double inTMin, double inTMax, double inStep);

private:
  typedef std::vector<State_T> X_Results_T_;
  struct StateAndTimeObserver;

  State_T x;
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
Results_T OdeintRunner<Rhs, Initializer>::run(double tMin, double tMax, double step) {
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
