#include "odeint_runner.h"

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
