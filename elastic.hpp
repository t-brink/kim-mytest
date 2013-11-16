/*
  Copyright (c) 2013 Tobias Brink

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include "core.hpp"

namespace mytest {
  /** Calculate energy using the Birch-Murnaghan equation of state. */
  double birch_murnaghan_energy(double V,
                                double E0, double V0,
                                double B0, double dB0_dp);
  /** Calculate energy using the Birch-Murnaghan equation of state. */
  inline
  double birch_murnaghan_energy(double V, BMParams p) {
    return birch_murnaghan_energy(V, p.E0, p.V0, p.B0, p.dB0_dp);
  }

  /** Calculate pressure using the Birch-Murnaghan equation of state. */
  double birch_murnaghan_pressure(double V, double V0,
                                  double B0, double dB0_dp);
  /** Calculate pressure using the Birch-Murnaghan equation of state. */
  inline
  double birch_murnaghan_pressure(double V, BMParams p) {
    return birch_murnaghan_pressure(V, p.V0, p.B0, p.dB0_dp);
  }

}


#endif /* ELASTIC_HPP */
