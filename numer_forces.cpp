/*
  Copyright (c) 2017 Tobias Brink

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

//#include "numer_forces.hpp"
#include "core.hpp"
#include "ndarray.hpp"

#include <limits>
#include <tuple>

using namespace std;
using namespace mytest;

double get_displ_energy(Compute& compute,
                        const unsigned i,
                        const unsigned dim,
                        const Vec3D<double>& orig_pos,
                        const double eps,
                        const std::map<std::string,int>& typemap) {
  Vec3D<double> new_pos(orig_pos);
  new_pos[dim] += eps;
  compute.set_position(i, new_pos, typemap);
  compute.update_neighbor_list();
  compute.compute();
  double E = compute.get_energy();
  compute.set_position(i, orig_pos, typemap);
  compute.update_neighbor_list();
  return E;
}


std::tuple<double, double, double>
df_dr(Compute& compute, unsigned i, unsigned dim, double eps,
      const std::map<std::string,int>& typemap) {
  static const unsigned ntab = 10;
  static const double con = 1.4;
  static const double con2 = con*con;
  static const double safe = 2.0;
  //
  Vec3D<double> orig_pos = compute.get_position(i);
  Array2D<double> a(ntab, ntab);
  double deriv = 0.0;
  double err = numeric_limits<double>::max();
  for (unsigned ii = 0; ii < ntab; ++ii) {
    double ep, en;
    // Displace and calculate energy (positive and negative direction).
    ep = get_displ_energy(compute, i, dim, orig_pos, +eps, typemap);
    en = get_displ_energy(compute, i, dim, orig_pos, -eps, typemap);
    a(0,ii) = (ep - en) / (2 * eps);
    double fac = con2;
    for (unsigned jj = 1; jj <= ii; ++jj) {
      a(jj, ii) = ( a(jj-1,ii) * fac - a(jj-1,ii-1) ) / (fac - 1);
      fac *= con2;
      double err_new = max( abs(a(jj,ii) - a(jj-1,ii)),
                            abs(a(jj,ii) - a(jj-1,ii-1)) );
      if (err_new < err) {
        deriv = a(jj,ii);
        err = err_new;
      }
    }
    if ( abs(a(ii,ii) - a(ii-1,ii-1)) > safe*err )
      break;
    eps /= con;
  }

  return make_tuple(deriv, eps, err);
}


// This guy computes forces numerically by differentiation and
// compares with the result of the model. Return value: maximum
// deviation.
double Compute::max_diff_force() {
  double max_diff = 0.0;
  // Calculate forces by model.
  Array2D<double> the_forces(box_->natoms, 3);
  this->compute();
  for (unsigned i = 0; i < box_->natoms; ++i) {
    for (unsigned dim = 0; dim < 3; ++dim) {
      the_forces(i, dim) = forces[3*i + dim];
    }
  }
  // Calculate forces by numerical derivatives.
  for (unsigned i = 0; i < box_->natoms; ++i) {
    for (unsigned dim = 0; dim < 3; ++dim) {
      double eps = 1e-6;
      double last_err = numeric_limits<double>::max();
      double err, deriv, last_deriv = 0.0;
      for (unsigned iteration = 0; iteration < 15; ++iteration) {
        auto retval = df_dr(*this, i, dim, eps, this->partcl_type_codes);
        deriv = get<0>(retval);
        eps = get<1>(retval);
        err = get<2>(retval);
        if (err > last_err) {
          deriv = last_deriv;
          err = last_err;
          break;
        }
        eps *= 10;
        last_deriv = deriv;
        last_err = err;
      }
      double num_force = -deriv;
      double diff = abs(the_forces(i, dim) - num_force);
      if (diff > max_diff)
        max_diff = diff;
    }
  }
  return max_diff;
}

