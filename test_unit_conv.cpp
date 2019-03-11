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

void print_it(const Compute& comp, double length_conv, double energy_conv) {
  printf("energy: %g\n", comp.get_energy() * energy_conv);
  printf("particle energies: %g %g %g\n",
         comp.get_energy(0) * energy_conv,
         comp.get_energy(1) * energy_conv,
         comp.get_energy(2) * energy_conv);
  printf("positions: %12g %12g %12g\n",
         comp.get_position(0)[0] * length_conv,
         comp.get_position(0)[1] * length_conv,
         comp.get_position(0)[2] * length_conv);
  printf("           %12g %12g %12g\n",
         comp.get_position(1)[0] * length_conv,
         comp.get_position(1)[1] * length_conv,
         comp.get_position(1)[2] * length_conv);
  printf("forces: %12g %12g %12g\n",
         comp.get_force(0,0) * energy_conv/length_conv,
         comp.get_force(0,1) * energy_conv/length_conv,
         comp.get_force(0,2) * energy_conv/length_conv);
  printf("        %12g %12g %12g\n",
         comp.get_force(1,0) * energy_conv/length_conv,
         comp.get_force(1,1) * energy_conv/length_conv,
         comp.get_force(1,2) * energy_conv/length_conv);
  printf("        %12g %12g %12g\n",
         comp.get_force(2,0) * energy_conv/length_conv,
         comp.get_force(2,1) * energy_conv/length_conv,
         comp.get_force(2,2) * energy_conv/length_conv);
  printf("virial:\n");
  printf("  %12g %12g %12g\n"
         "               %12g %12g\n"
         "                            %12g\n",
         comp.get_virial().xx * energy_conv,
         comp.get_virial().xy * energy_conv,
         comp.get_virial().xz * energy_conv,
         comp.get_virial().yy * energy_conv,
         comp.get_virial().yz * energy_conv,
         comp.get_virial().zz * energy_conv);
  printf("virial from dE/dr:\n");
  printf("  %12g %12g %12g\n"
         "               %12g %12g\n"
         "                            %12g\n",
         comp.get_virial_from_dEdr().xx * energy_conv,
         comp.get_virial_from_dEdr().xy * energy_conv,
         comp.get_virial_from_dEdr().xz * energy_conv,
         comp.get_virial_from_dEdr().yy * energy_conv,
         comp.get_virial_from_dEdr().yz * energy_conv,
         comp.get_virial_from_dEdr().zz * energy_conv);
}

int main() {
  // TODO                                         
  // * do not print for the user to compare, but compare like in
  //   test_switch_params automatically
  // * more potentials


  // Create the simulation cell. ///////////////////////////////////////
  random_device rd;
  mt19937 rng(rd());
  vector<string> atomtypes = { "Si", "C" };
  const auto random_box = Box::random_box(10.0, 10.0, 10.0,
                                          true, true, false,
                                          1.75, atomtypes, "boxname", rng);

  // Now Erhart/Albe ///////////////////////////////////////////////////
  Compute erhart2005(make_unique<Box>(*random_box),
                     "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_003");
  erhart2005.compute();

  print_it(erhart2005, 1.0, 1.0);


  // Convert lengths to nm.
  Compute erhart2005_nm(make_unique<Box>(*random_box),
                        "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_003",
                        KIM::LENGTH_UNIT::nm,
                        KIM::ENERGY_UNIT::eV);
  erhart2005_nm.scale_box(0.1);
  erhart2005_nm.update_neighbor_list();
  erhart2005_nm.compute();

  print_it(erhart2005_nm, 10.0, 1.0);


  // Convert energy to J.
  Compute erhart2005_J(make_unique<Box>(*random_box),
                       "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_003",
                       KIM::LENGTH_UNIT::A,
                       KIM::ENERGY_UNIT::J);
  erhart2005_J.compute();

  print_it(erhart2005_J, 1.0, 6.241508e18);

}
