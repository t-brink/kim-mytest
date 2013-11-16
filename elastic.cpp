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

#include "elastic.hpp"

#include <nlopt.hpp>

#include <utility>

using namespace std;
using namespace mytest;

static
double birch_murnaghan_energy(double V,
                              double E0, double V0,
                              double B0, double dB0_dp) {
  const double V_ = pow(V0 / V, 2.0/3.0);
  return E0 + 9*V0*B0/16 * (pow3(V_ - 1) * dB0_dp
                            + pow2(V_ - 1) * (6 - 4*V_));
}

static
double obj_func_birch_murnaghan_energy(const vector<double>& x,
                                       vector<double>& grad,
                                       void* f_data) {
  // Reference volumes and energies.
  auto& VE = *static_cast<pair<vector<double>&,vector<double>&>*>(f_data);
  // Gradient.
  if (!grad.empty())
    throw runtime_error("Gradient not supported.");
  double obj_func = 0.0;
  for(unsigned i = 0; i != VE.first.size(); ++i) {
    const double V = VE.first[i];
    const double ref_E = VE.second[i];
    obj_func += abs(ref_E - birch_murnaghan_energy(V, x[0], x[1], x[2], x[3]));
  }
  return obj_func;
}

BMParams mytest::bulk_modulus_energy(Compute& compute,
                                     bool c_to_a, bool b_to_a,
                                     bool positions,
                                     bool angle_ab,
                                     bool angle_ac,
                                     bool angle_bc) {
  if (c_to_a || b_to_a)
    throw runtime_error("not implemented");
  if (positions)
    throw runtime_error("not implemented");
  if (angle_ab || angle_ac || angle_bc)
    throw runtime_error("not implemented");
  // Set up 11 boxes with different volumes.
  const double epsilon = 0.01; // Total maximum strain is 5%.
  vector<double> volumes;
  vector<unique_ptr<Box>> boxes;
  for (int i = -5; i <= 5; ++i) {
    unique_ptr<Box> p = compute.copy_box();
    p->scale(i*epsilon + 1);
    volumes.push_back(p->calc_volume());
    boxes.push_back(move(p));
  }
  // Switch through boxes, recording energy.
  vector<double> energies;
  // Last entry of boxes now contains the original box.  As we always
  // calculate energy for the (i-1)th element and then put the ith
  // element into compute, we will end up with the original box in the
  // Compute object without calculating energy for it once.  The last
  // element of boxes will be a null pointer after the loop.
  boxes.push_back(compute.change_box(move(boxes[0])));
  for (unsigned i = 1; i != boxes.size(); ++i) {
    compute.compute();
    energies.push_back(compute.get_energy());
    boxes[i-1] = compute.change_box(move(boxes[i]));
  }

  // Now fit Birch-Murnaghan.
  double obj_val;
  vector<double> parameters = { // Inital guess.
    energies[5], volumes[5], 1.0, 1.0
  };
  pair<vector<double>&,vector<double>&> refdata(volumes, energies);
  nlopt::opt optimizer(nlopt::LN_SBPLX, 4);
  optimizer.set_min_objective(obj_func_birch_murnaghan_energy, &refdata);
  optimizer.set_ftol_abs(1e-6); // TODO: too much/too little?
  optimizer.optimize(parameters, obj_val);

  /*
  for (int i = 0; i != volumes.size(); ++i)
    cout << volumes[i] << "   "
         << energies[i] << "   "
         << birch_murnaghan_energy(volumes[i],
                                   parameters[0],
                                   parameters[1],
                                   parameters[2],
                                   parameters[3]) << endl;
  */

  return { parameters[0], parameters[1], parameters[2], parameters[3] };
}
