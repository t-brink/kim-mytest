/*
  Copyright (c) 2018 Tobias Brink

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

#include <iostream>
#include <cstdio>

#include "utils.hpp"
#include "core.hpp"
#include "elastic.hpp"

using namespace std;
using namespace mytest;

struct AllResults {
  double energy;
  vector<double> particle_energy;
  vector<double> force;
  Voigt6<double> virial;
  vector<Voigt6<double>> particle_virial;

  bool approx_equal(const AllResults&, double);
};

bool AllResults::approx_equal(const AllResults& other, double eps = 1e-6) {
  if (abs(energy - other.energy) > eps)
    return false;
  if (particle_energy.size() != other.particle_energy.size())
    return false;
  for (unsigned i = 0; i < particle_energy.size(); ++i) {
    if (abs(particle_energy[i] - other.particle_energy[i]) > eps)
      return false;
  }
  if (force.size() != other.force.size())
    return false;
  for (unsigned i = 0; i < force.size(); ++i) {
    if (abs(force[i] - other.force[i]) > eps)
      return false;
  }
  for (unsigned dim = 0; dim < 6; ++dim) {
    if (abs(virial(dim) - other.virial(dim)) > eps)
      return false;
  }
  if (particle_virial.size() != other.particle_virial.size())
    return false;
  for (unsigned i = 0; i < particle_virial.size(); ++i) {
    for (unsigned dim = 0; dim < 6; ++dim) {
      if (abs(particle_virial[i](dim) - other.particle_virial[i](dim)) > eps)
        return false;
    }
  }
  return true;
}

AllResults compute(Compute& comp) {
  comp.compute();
  unsigned natoms = comp.get_natoms();
  AllResults res;
  res.energy = comp.get_energy();
  for (unsigned dim = 0; dim < 6; ++dim) {
    res.virial(dim) = comp.get_virial()(dim);
  }
  for (unsigned i = 0; i < natoms; ++i) {
    res.particle_energy.push_back(comp.get_energy(i));
    for (unsigned dim = 0; dim < 3; ++dim) {
      res.force.push_back(comp.get_force(i, dim));
    }
    res.particle_virial.push_back(comp.get_virial(i));
  }
  return move(res);
}

void do_the_shit(unique_ptr<Box> box, const string& modelname,
                 const AllResults& res_tersoff1988,
                 const AllResults& res_tersoff1989,
                 const AllResults& res_erhart2005) {
  Compute comp(move(box), modelname);
  cout << "model: " << modelname << endl;
  /*
  const vector<unsigned> idx2{ 0, 0 };
  cout << "    A  = " << comp.get_parameter_double("PARAM_FREE_A", idx2) << endl;
  cout << "    B  = " << comp.get_parameter_double("PARAM_FREE_B", idx2) << endl;
  cout << "    λ₁ = " << comp.get_parameter_double("PARAM_FREE_lambda1", idx2) << endl;
  cout << "    λ₂ = " << comp.get_parameter_double("PARAM_FREE_lambda2", idx2) << endl;
  cout << "    β  = " << comp.get_parameter_double("PARAM_FREE_beta", idx2) << endl;
  cout << "    n  = " << comp.get_parameter_double("PARAM_FREE_n", idx2) << endl;
  const vector<unsigned> idx3{ 0, 0, 0 };
  cout << "    λ₃ = " << comp.get_parameter_double("PARAM_FREE_lambda3", idx3) << endl;
  cout << "    m  = " << comp.get_parameter_int("PARAM_FREE_m", idx3) << endl;
  cout << "    γ  = " << comp.get_parameter_double("PARAM_FREE_gamma", idx3) << endl;
  cout << "    c  = " << comp.get_parameter_double("PARAM_FREE_c", idx3) << endl;
  cout << "    d  = " << comp.get_parameter_double("PARAM_FREE_d", idx3) << endl;
  cout << "    h  = " << comp.get_parameter_double("PARAM_FREE_h", idx3) << endl;
  cout << "    Rc = " << comp.get_parameter_double("PARAM_FREE_Rc", idx3) << endl;
  cout << "    Dc = " << comp.get_parameter_double("PARAM_FREE_Dc", idx3) << endl;
  */
  comp.compute();
  cout << "E = " << comp.get_energy() << " eV" << endl;

  // Set to Tersoff PRB37 1988 (all types Si)
  for (unsigned i = 0; i < comp.get_number_of_particle_types(); ++i) {
  for (unsigned j = 0; j < comp.get_number_of_particle_types(); ++j) {
  for (unsigned k = 0; k < comp.get_number_of_particle_types(); ++k) {
  const vector<unsigned> idx2{ i, j };
  const vector<unsigned> idx3{ i, j, k };
  comp.set_parameter("PARAM_FREE_A", idx2, 3264.7, false);
  comp.set_parameter("PARAM_FREE_B", idx2, 95.373, false);
  comp.set_parameter("PARAM_FREE_lambda1", idx2, 3.2394, false);
  comp.set_parameter("PARAM_FREE_lambda2", idx2, 1.3258, false);
  comp.set_parameter("PARAM_FREE_beta", idx2, 0.33675, false);
  comp.set_parameter("PARAM_FREE_n", idx2, 22.956, false);
  comp.set_parameter("PARAM_FREE_lambda3", idx3, 1.3258, false);
  comp.set_parameter("PARAM_FREE_m", idx3, 3, false);
  comp.set_parameter("PARAM_FREE_gamma", idx3, 1.0, false);
  comp.set_parameter("PARAM_FREE_c", idx3, 4.8381, false);
  comp.set_parameter("PARAM_FREE_d", idx3, 2.0417, false);
  comp.set_parameter("PARAM_FREE_h", idx3, 0.0, false);
  comp.set_parameter("PARAM_FREE_Rc", idx3, 3.0, false);
  comp.set_parameter("PARAM_FREE_Dc", idx3, 0.2, true);
  }}}
  comp.update_neighbor_list();
  AllResults new_tersoff1988 = compute(comp);
  cout << "Tersoff PRB37 1988: ";
  if (new_tersoff1988.approx_equal(res_tersoff1988)) {
    cout << "OK";
  } else {
    cout << "mismatch!";
  }
  cout << endl;

  // Set to Tersoff PRB39 1989 (all types Si)
  for (unsigned i = 0; i < comp.get_number_of_particle_types(); ++i) {
  for (unsigned j = 0; j < comp.get_number_of_particle_types(); ++j) {
  for (unsigned k = 0; k < comp.get_number_of_particle_types(); ++k) {
  const vector<unsigned> idx2{ i, j };
  const vector<unsigned> idx3{ i, j, k };
  comp.set_parameter("PARAM_FREE_A", idx2, 1830.8, false);
  comp.set_parameter("PARAM_FREE_B", idx2, 471.18, false);
  comp.set_parameter("PARAM_FREE_lambda1", idx2, 2.4799, false);
  comp.set_parameter("PARAM_FREE_lambda2", idx2, 1.7322, false);
  comp.set_parameter("PARAM_FREE_beta", idx2, 1.1e-06, false);
  comp.set_parameter("PARAM_FREE_n", idx2, 0.78734, false);
  comp.set_parameter("PARAM_FREE_lambda3", idx3, 0.0, false);
  comp.set_parameter("PARAM_FREE_m", idx3, 3, false);
  comp.set_parameter("PARAM_FREE_gamma", idx3, 1.0, false);
  comp.set_parameter("PARAM_FREE_c", idx3, 100390.0, false);
  comp.set_parameter("PARAM_FREE_d", idx3, 16.217, false);
  comp.set_parameter("PARAM_FREE_h", idx3, -0.59825, false);
  comp.set_parameter("PARAM_FREE_Rc", idx3, 2.85, false);
  comp.set_parameter("PARAM_FREE_Dc", idx3, 0.15, true);
  }}}
  comp.update_neighbor_list();
  AllResults new_tersoff1989 = compute(comp);
  cout << "Tersoff PRB39 1989: ";
  if (new_tersoff1989.approx_equal(res_tersoff1989)) {
    cout << "OK";
  } else {
    cout << "mismatch!";
  }
  cout << endl;

  // Set to Erhart/Albe (all types Si)
  for (unsigned i = 0; i < comp.get_number_of_particle_types(); ++i) {
  for (unsigned j = 0; j < comp.get_number_of_particle_types(); ++j) {
  for (unsigned k = 0; k < comp.get_number_of_particle_types(); ++k) {
  const vector<unsigned> idx2{ i, j };
  const vector<unsigned> idx3{ i, j, k };
  comp.set_parameter("PARAM_FREE_A", idx2, 2145.7128, false);
  comp.set_parameter("PARAM_FREE_B", idx2, 219.521624, false);
  comp.set_parameter("PARAM_FREE_lambda1", idx2, 2.83318929, false);
  comp.set_parameter("PARAM_FREE_lambda2", idx2, 1.53810493, false);
  comp.set_parameter("PARAM_FREE_beta", idx2, 1.0, false);
  comp.set_parameter("PARAM_FREE_n", idx2, 1.0, false);
  comp.set_parameter("PARAM_FREE_lambda3", idx3, 0.0, false);
  comp.set_parameter("PARAM_FREE_m", idx3, 1, false);
  comp.set_parameter("PARAM_FREE_gamma", idx3, 0.114354, false);
  comp.set_parameter("PARAM_FREE_c", idx3, 2.00494, false);
  comp.set_parameter("PARAM_FREE_d", idx3, 0.81472, false);
  comp.set_parameter("PARAM_FREE_h", idx3, -0.259, false);
  comp.set_parameter("PARAM_FREE_Rc", idx3, 2.82, false);
  comp.set_parameter("PARAM_FREE_Dc", idx3, 0.14, true);
  }}}
  comp.update_neighbor_list();
  AllResults new_erhart2005 = compute(comp);
  cout << "Erhart/Albe 2005  : ";
  if (new_erhart2005.approx_equal(res_erhart2005)) {
    cout << "OK";
  } else {
    cout << "mismatch!  " << res_erhart2005.energy << " <-> " << new_erhart2005.energy;
  }
  cout << endl;
}

int main() {
  // RNG
  random_device rd;
  mt19937 rng(rd());
  // Create the simulation cell.
  const auto rand_box = Box::random_box(20.0, 20.0, 20.0, true, false, false, 2.1,
                                        "Si", KIM_neigh_pure_f, "mybox", rng);
  // Set up Tersoff PRB37 1988 ////////////////////////////////////////////////////////////////////
  Compute tersoff1988(make_unique<Box>(*rand_box),
                      "Tersoff_LAMMPS_Tersoff_PRB37_1988_Si__MO_245095684871_000");
  AllResults res_tersoff1988 = compute(tersoff1988);

  // Now Tersoff PRB39 1989 ///////////////////////////////////////////////////////////////////////
  Compute tersoff1989(make_unique<Box>(*rand_box),
                      "Tersoff_LAMMPS_Tersoff_PRB39_1989_CSi__MO_171585019474_000");
  AllResults res_tersoff1989 = compute(tersoff1989);

  // Now Erhart/Albe //////////////////////////////////////////////////////////////////////////////
  Compute erhart2005(make_unique<Box>(*rand_box),
                     "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002");
  AllResults res_erhart2005 = compute(erhart2005);

  // Test it ///////////////////////////////////////////////////////////////////////////////////////
  do_the_shit(make_unique<Box>(*rand_box),
              "Tersoff_LAMMPS_Tersoff_PRB37_1988_Si__MO_245095684871_000",
              res_tersoff1988, res_tersoff1989, res_erhart2005);

  cout << endl;

  do_the_shit(make_unique<Box>(*rand_box),
              "Tersoff_LAMMPS_Tersoff_PRB39_1989_CSi__MO_171585019474_000",
              res_tersoff1988, res_tersoff1989, res_erhart2005);

  cout << endl;

  do_the_shit(make_unique<Box>(*rand_box),
              "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002",
              res_tersoff1988, res_tersoff1989, res_erhart2005);
}
