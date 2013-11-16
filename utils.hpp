#ifndef MYTEST_UTILS_HPP
#define MYTEST_UTILS_HPP

#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

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

  /** Return y-intercept and slope obtained by linear regression. */
  inline
  std::pair<double,double> linear_leastsq(const std::vector<double>& x,
                                          const std::vector<double>& y){
    if (x.size() != y.size())
      throw std::runtime_error("Input vectors must have same length");
    const double av_x = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    const double av_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double slope_numerator = 0.0, slope_denominator = 0.0;
    for (unsigned i = 0; i != x.size(); ++i) {
      const double dx = x[i] - av_x;
      slope_numerator += dx * (y[i] - av_y);
      slope_denominator += dx*dx;
    }
    const double slope = slope_numerator / slope_denominator;
    const double intercept = av_y - slope*av_x;
    return { intercept, slope };
  }
}


#endif /* MYTEST_UTILS_HPP */
