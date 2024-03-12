#ifndef ODEINT_RUNNER_H
#define ODEINT_RUNNER_H

#include <boost/numeric/odeint.hpp>

#include <vector>

// For State_T x, x[0] represents f(t) and x[1] represents f'(t).
typedef std::array<double, 2> State_T;

typedef std::vector<double> ResultSeq_T;
typedef ResultSeq_T TimeSeq_T;

typedef std::pair<ResultSeq_T, TimeSeq_T> Results_T;

namespace detail
{

struct StateAndTimeObserver
{
  std::vector<State_T> &m_states;
  std::vector<double> &m_times;

  StateAndTimeObserver(std::vector<State_T> &states, std::vector<double> &times)
      : m_states(states), m_times(times)
  {}

  void operator()(const State_T &x, double t)
  {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

} // namespace detail

template <typename T>
concept RHS = requires(T a, const State_T &x, State_T &dxdt, const double t) { a(x, dxdt, t); };

template <typename T>
concept Initializer = requires(T a) {
  {
    a.init()
  } -> std::same_as<State_T>;
};

template <RHS Rhs, Initializer Init>
class OdeintRunner
{
public:
  OdeintRunner() { x = Init{}.init(); }

  Results_T run(double inTMin, double inTMax, double inStep);

private:
  typedef std::vector<State_T> X_Results_T_;

  State_T x{};
};

template <RHS Rhs, Initializer Init>
Results_T OdeintRunner<Rhs, Init>::run(double tMin, double tMax, double step)
{
  namespace boostode = boost::numeric::odeint;

  X_Results_T_ x_vec;
  TimeSeq_T times;

  // With constant stepper and integrator, observer is called at regular intervals.
  boostode::runge_kutta4<State_T> stepper;
  boostode::integrate_const(stepper, Rhs{}, x, tMin, tMax, step, detail::StateAndTimeObserver(x_vec, times));

  ResultSeq_T resultsFOnly{};
  std::transform(begin(x_vec), end(x_vec), std::back_inserter(resultsFOnly),
                 [](const State_T &x) { return x[0]; });

  return std::make_pair(std::move(resultsFOnly), std::move(times));
}

#endif
