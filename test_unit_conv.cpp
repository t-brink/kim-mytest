/*
  Copyright (c) 2019 Tobias Brink

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
#include <random>

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
  Voigt6<double> virial_from_dEdr;
  vector<Voigt6<double>> particle_virial;

  bool approx_equal(const AllResults&, double);
};

bool AllResults::approx_equal(const AllResults& other, double eps = 1e-6) {
  if (abs(energy - other.energy) > eps) {
    cout << "ENERGY ";
    return false;
  }
  if (particle_energy.size() != other.particle_energy.size()) {
    cout << "ATOM ENERGY SIZE ";
    return false;
  }
  for (unsigned i = 0; i < particle_energy.size(); ++i) {
    if (abs(particle_energy[i] - other.particle_energy[i]) > eps) {
      cout << "ATOM ENERGY ";
      return false;
    }
  }
  if (force.size() != other.force.size()) {
    cout << "FORCE SIZE ";
    return false;
  }
  for (unsigned i = 0; i < force.size(); ++i) {
    if (abs(force[i] - other.force[i]) > eps) {
      cout << "FORCE " << force[i] << " " << other.force[i] << " ";
      return false;
    }
  }
  for (unsigned dim = 0; dim < 6; ++dim) {
    if (abs(virial(dim) - other.virial(dim)) > 10*eps) {
      cout << "VIRIAL ";
      return false;
    }
  }
  for (unsigned dim = 0; dim < 6; ++dim) {
    if (abs(virial_from_dEdr(dim) - other.virial_from_dEdr(dim)) > 10*eps) {
      cout << "VIRIAL FROM dE/dr ";
      return false;
    }
  }
  if (particle_virial.size() != other.particle_virial.size()) {
    cout << "ATOM VIRIAL SIZE ";
    return false;
  }
  for (unsigned i = 0; i < particle_virial.size(); ++i) {
    for (unsigned dim = 0; dim < 6; ++dim) {
      if (abs(particle_virial[i](dim) - other.particle_virial[i](dim)) > 10*eps) {
        cout << "ATOM VIRIAL ";
        return false;
      }
    }
  }
  return true;
}

AllResults compute(Compute& comp, double length_conv, double energy_conv) {
  // If given the correct args, convert everything to eV/Å.
  comp.compute();
  unsigned natoms = comp.get_natoms();
  AllResults res;
  res.energy = comp.get_energy() * energy_conv;
  for (unsigned dim = 0; dim < 6; ++dim) {
    res.virial(dim) = comp.get_virial()(dim) * energy_conv;
    res.virial_from_dEdr(dim) = comp.get_virial_from_dEdr()(dim) * energy_conv;
  }
  for (unsigned i = 0; i < natoms; ++i) {
    res.particle_energy.push_back(comp.get_energy(i) * energy_conv);
    for (unsigned dim = 0; dim < 3; ++dim) {
      res.force.push_back(comp.get_force(i, dim) * energy_conv/length_conv);
    }
    res.particle_virial.push_back(comp.get_virial(i) * energy_conv);
  }
  return move(res);
}



void do_it(unique_ptr<Box> box, const string& model, const string& modelname) {
  Compute comp(make_unique<Box>(*box), model);
  AllResults res = compute(comp, 1.0, 1.0);


  // Convert lengths to nm.
  Compute comp_nm(make_unique<Box>(*box), model,
                  KIM::LENGTH_UNIT::nm,
                  KIM::ENERGY_UNIT::eV);
  comp_nm.scale_box(0.1);
  comp_nm.update_neighbor_list();
  AllResults res_nm = compute(comp_nm, 10.0, 1.0);


  // Convert energy to J.
  Compute comp_J(make_unique<Box>(*box), model,
                 KIM::LENGTH_UNIT::A,
                 KIM::ENERGY_UNIT::J);
  AllResults res_J = compute(comp_J, 1.0, 6.241508e18);


  // Convert length to nm and energy to J.
  Compute comp_nm_J(make_unique<Box>(*box), model,
                    KIM::LENGTH_UNIT::nm,
                    KIM::ENERGY_UNIT::J);
  comp_nm_J.scale_box(0.1);
  comp_nm_J.update_neighbor_list();
  AllResults res_nm_J = compute(comp_nm_J, 10.0, 6.241508e18);

  // Compare.
  cout << modelname << " (eV/Å) : "
       << res.energy
       << " eV (reference)"
       << endl;


  cout << modelname << " (eV/nm): "
       << res_nm.energy
       << " eV ";
  if (res.approx_equal(res_nm, 1e-6)) {
    cout << "OK";
  } else {
    cout << "mismatch!!!";
  }
  cout << endl;


  cout << modelname << " (J/Å)  : "
       << res_J.energy
       << " eV ";
  if (res.approx_equal(res_J, 1e-4)) {
    cout << "OK";
  } else {
    cout << "mismatch!!!";
  }
  cout << endl;


  cout << modelname << " (J/nm) : "
       << res_nm_J.energy
       << " eV ";
  if (res.approx_equal(res_nm_J, 1e-4)) {
    cout << "OK";
  } else {
    cout << "mismatch!!!";
  }
  cout << endl;

}



int main() {
  random_device rd;
  mt19937 rng(rd());

  // Create simulation cells. //////////////////////////////////////////
  vector<string> atomtypes_Si = { "Si" };
  const auto Si_box = Box::random_box(10.0, 10.0, 10.0,
                                      true, true, false,
                                      2.1, atomtypes_Si, "boxname", rng);

  vector<string> atomtypes_SiC = { "Si", "C" };
  const auto SiC_box = Box::random_box(10.0, 10.0, 10.0,
                                       true, true, false,
                                       1.75, atomtypes_SiC, "boxname", rng);

  vector<string> atomtypes_PtC = { "Pt", "C" };
  const auto PtC_box = Box::random_box(10.0, 10.0, 10.0,
                                       true, true, false,
                                       2.5, atomtypes_PtC, "boxname", rng);

  // Tersoff PRB39 1989 ////////////////////////////////////////////////
  do_it(make_unique<Box>(*Si_box),
        "Tersoff_LAMMPS_Tersoff_PRB37_1988_Si__MO_245095684871_001",
        "Tersoff PRB37 1988");

  cout << endl;

  // Tersoff PRB39 1989 ////////////////////////////////////////////////
  do_it(make_unique<Box>(*SiC_box),
        "Tersoff_LAMMPS_Tersoff_PRB39_1989_CSi__MO_171585019474_001",
        "Tersoff PRB39 1989");

  cout << endl;

  // Now Erhart/Albe ///////////////////////////////////////////////////
  do_it(make_unique<Box>(*SiC_box),
        "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_003",
        "Erhart/Albe 2005");

  cout << endl;

  // Now Erhart/Albe ///////////////////////////////////////////////////
  do_it(make_unique<Box>(*PtC_box),
        "Tersoff_LAMMPS_Albe_Nordlund_Averback_PtC__MO_500121566391_001",
        "Albe et al. 2002");
}
