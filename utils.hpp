#ifndef MYTEST_UTILS_HPP
#define MYTEST_UTILS_HPP

#include <cmath>
#include <memory>

namespace mytest {
  /** Modulo with the same semantics as python (or maths). */
  inline
  double pmod(double x, double N) {
    return (x < 0) ? fmod(fmod(x, N) + N, N) : fmod(x, N);
  }

  /** Modulo with the same semantics as python (or maths).

      Will only actually do the modulo if the last argument is @c
      true.
  */
  inline
  double pmodif(double x, double N, bool do_mod) {
    if (do_mod)
      return (x < 0) ? fmod(fmod(x, N) + N, N) : fmod(x, N);
    else
      return x;
  }

  /** Utility function that should be in C++14 but not yet in C++11. */
  template<typename T, typename ...Args>
  std::unique_ptr<T> make_unique(Args&& ...args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
}


#endif /* MYTEST_UTILS_HPP */
