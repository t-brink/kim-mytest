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
#include "utils.hpp"

#include <nlopt.hpp>

#include <utility>

using namespace std;
using namespace mytest;

double mytest::birch_murnaghan_energy(double V,
                                      double E0, double V0,
                                      double B0, double dB0_dp) {
  const double V_ = pow(V0 / V, 2.0/3.0);
  return E0 + 9*V0*B0/16 * (pow3(V_ - 1) * dB0_dp
                            + pow2(V_ - 1) * (6 - 4*V_));
}

double mytest::birch_murnaghan_pressure(double V, double V0,
                                        double B0, double dB0_dp) {
  const double V7_3 = pow(V0 / V, 7.0/3.0);
  const double V5_3 = pow(V0 / V, 5.0/3.0);
  const double V2_3 = pow(V0 / V, 2.0/3.0);
  return 3*B0/2 * (V7_3 - V5_3) * (1 + 0.75 * (dB0_dp - 4) * (V2_3 - 1));
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

BMParams Compute::bulk_modulus_energy(vector<double>& volumes,
                                      vector<double>& energies,
                                      double max_strain,
                                      bool c_to_a, bool b_to_a,
                                      bool positions,
                                      bool angle_ab,
                                      bool angle_ac,
                                      bool angle_bc) {
  // Check if energy calculation is supported.
  if (!has_energy)
    throw runtime_error("The model does not provide energy calculation "
                        "which is needed to calculate bulk modulus.");
  //
  static const int n_boxes = 5; // Number of boxes in positive and
                                // negative strain direction (i.e. the
                                // real number of boxes is 2*n_boxes + 1).
  if (c_to_a || b_to_a)
    throw runtime_error("not implemented");
  if (angle_ab || angle_ac || angle_bc)
    throw runtime_error("not implemented");
  // Set up 11 boxes with different volumes.
  const double epsilon = max_strain / n_boxes;
  volumes.clear();
  vector<unique_ptr<Box>> boxes;
  for (int i = -5; i <= 5; ++i) {
    unique_ptr<Box> p = copy_box();
    p->scale(i*epsilon + 1);
    volumes.push_back(p->calc_volume());
    boxes.push_back(move(p));
  }
  // Switch through boxes, recording energy.
  energies.clear();
  // Last entry of boxes now contains the original box.  As we always
  // calculate energy for the (i-1)th element and then put the ith
  // element into compute, we will end up with the original box in the
  // Compute object without calculating energy for it once.  The last
  // element of boxes will be a null pointer after the loop.
  boxes.push_back(change_box(move(boxes[0])));
  for (unsigned i = 1; i != boxes.size(); ++i) {
    if (positions)
      optimize_positions(0.001, 10000);
    compute();
    energies.push_back(get_energy());
    boxes[i-1] = change_box(move(boxes[i]));
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

  return { parameters[0], parameters[1], parameters[2], parameters[3] };
}


static
double obj_func_birch_murnaghan_pressure(const vector<double>& x,
                                         vector<double>& grad,
                                         void* f_data) {
  // Reference volumes and pressures.
  auto& Vp = *static_cast<pair<vector<double>&,vector<double>&>*>(f_data);
  // Gradient.
  if (!grad.empty())
    throw runtime_error("Gradient not supported.");
  double obj_func = 0.0;
  for(unsigned i = 0; i != Vp.first.size(); ++i) {
    const double V = Vp.first[i];
    const double ref_p = Vp.second[i];
    obj_func += abs(ref_p - birch_murnaghan_pressure(V, x[0], x[1], x[2]));
  }
  return obj_func;
}

BMParams Compute::bulk_modulus_pressure(vector<double>& volumes,
                                        vector<double>& pressures,
                                        double max_strain,
                                        bool c_to_a, bool b_to_a,
                                        bool positions,
                                        bool angle_ab,
                                        bool angle_ac,
                                        bool angle_bc) {
  // Check if virial calculation is supported.
  if (!has_virial)
    throw runtime_error("The model does not provide virial calculation "
                        "which is needed to calculate bulk modulus.");
  //
  static const int n_boxes = 5; // Number of boxes in positive and
                                // negative strain direction (i.e. the
                                // real number of boxes is 2*n_boxes + 1).
  if (c_to_a || b_to_a)
    throw runtime_error("not implemented");
  if (angle_ab || angle_ac || angle_bc)
    throw runtime_error("not implemented");
  // Set up 11 boxes with different volumes.
  volumes.clear();
  const double epsilon = max_strain / n_boxes;
  vector<unique_ptr<Box>> boxes;
  for (int i = -5; i <= 5; ++i) {
    unique_ptr<Box> p = copy_box();
    p->scale(i*epsilon + 1);
    volumes.push_back(p->calc_volume());
    boxes.push_back(move(p));
  }
  // Switch through boxes, recording pressure.
  pressures.clear();
  // Last entry of boxes now contains the original box.  As we always
  // calculate pressure for the (i-1)th element and then put the ith
  // element into compute, we will end up with the original box in the
  // Compute object without calculating pressure for it once.  The last
  // element of boxes will be a null pointer after the loop.
  boxes.push_back(change_box(move(boxes[0])));
  for (unsigned i = 1; i != boxes.size(); ++i) {
    if (positions)
      optimize_positions(0.001, 10000);
    compute();
    // Assume hydrostatic, take average.
    const Voigt6<double> press = get_virial();
    pressures.push_back(-(press.xx + press.yy + press.zz) / (3 * volumes[i-1]));
    boxes[i-1] = change_box(move(boxes[i]));
  }

  // Now fit Birch-Murnaghan.
  double obj_val;
  vector<double> parameters = { // Inital guess.
    volumes[5], 1.0, 1.0
  };
  pair<vector<double>&,vector<double>&> refdata(volumes, pressures);
  nlopt::opt optimizer(nlopt::LN_SBPLX, 3);
  optimizer.set_min_objective(obj_func_birch_murnaghan_pressure, &refdata);
  optimizer.set_ftol_abs(1e-6); // TODO: too much/too little?
  optimizer.optimize(parameters, obj_val);

  return { 0.0, parameters[0], parameters[1], parameters[2] };
}

double Compute::elastic_constant(unsigned i, unsigned j,
                                 vector<double>& strain,
                                 vector<double>& stress,
                                 double max_strain,
                                 bool positions) {
  // Check if virial calculation is supported.
  if (!has_virial)
    throw runtime_error("The model does not provide virial calculation "
                        "which is needed to calculate elastic constants.");
  //
  static const int n_boxes = 5; // Number of boxes in positive and
                                // negative strain direction (i.e. the
                                // real number of boxes is 2*n_boxes + 1).
  if (i < 1 || i > 6)
    throw runtime_error("i must be in [1..6]");
  if (j < 1 || j > 6)
    throw runtime_error("j must be in [1..6]");
  const unsigned ii = i - 1;
  const unsigned jj = j - 1;
  // Set up 11 boxes with different strains.
  strain.clear();
  vector<double> volumes;
  const double epsilon = max_strain / n_boxes;
  vector<unique_ptr<Box>> boxes;
  for (int i = -5; i <= 5; ++i) {
    unique_ptr<Box> p = copy_box();
    const double eps = i*epsilon;
    Voigt6<double> defmatrix(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    defmatrix(ii) += eps;
    p->deform(defmatrix);
    strain.push_back(eps);
    volumes.push_back(p->calc_volume());
    boxes.push_back(move(p));
  }
  // Switch through boxes, recording stress.
  stress.clear();
  // Last entry of boxes now contains the original box.  As we always
  // calculate pressure for the (i-1)th element and then put the ith
  // element into compute, we will end up with the original box in the
  // Compute object without calculating pressure for it once.  The last
  // element of boxes will be a null pointer after the loop.
  boxes.push_back(change_box(move(boxes[0])));
  for (unsigned i = 1; i != boxes.size(); ++i) {
    if (positions)
      optimize_positions(0.001, 10000);
    compute();
    const Voigt6<double> press = get_virial();
    stress.push_back(press(jj) / volumes[i-1]);
    boxes[i-1] = change_box(move(boxes[i]));
  }

  // Now fit Hooke.
  auto res = linear_leastsq(strain, stress);
  const double eps0 = res.first;
  const double c = (j <= 3) ? res.second : res.second/2;
  if (abs(eps0) > 0.01)
    cerr << "WARNING: stress at strain = 0 is not zero, "
         << "box may not be well optimized." << endl;
  return c;
}
