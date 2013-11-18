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

#include <iostream>
#include <cstdio>

#include "utils.hpp"
#include "core.hpp"
#include "elastic.hpp"
#include "dsl.hpp"

using namespace std;
using namespace mytest;

int main() {
  map< string,unique_ptr<Box> > boxes;
  map<string,Compute> computes;
  string input;
  while (getline(cin, input))
    parse(input, boxes, computes);
  return 0;

  /*
  const vector<int> types1{ 0 };
  Compute comp1(make_unique<Box>("fcc", 3.940, true, 3, 3, 3,
                                 true, true, true,
                                 types1, KIM_neigh_rvec_f, "box1"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  const vector<int> types2{ 0, 1 };
  Compute comp2(make_unique<Box>("B3", 4.359, true, 3, 3, 3,
                                 true, true, true,
                                 types2, KIM_neigh_rvec_f, "box2"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");

  comp1.compute();
  comp2.compute();
  cout << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom    "
       << "comp2 = " << comp2.get_energy_per_atom() << " eV/atom" << endl;

  comp1.switch_boxes(comp2);

  comp1.compute();
  comp2.compute();
  cout << "comp2 = " << comp2.get_energy_per_atom() << " eV/atom    "
       << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom" << endl;

  Compute comp3(make_unique<Box>("B1", 4.244, true, 3, 3, 3,
                                 true, true, true,
                                 types2, KIM_neigh_rvec_f, "box"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  Compute comp4(make_unique<Box>("B1", 4.244, false, 3, 3, 3,
                                 true, true, true,
                                 types2, KIM_neigh_rvec_f, "box"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  Compute comp5(make_unique<Box>("B2", 2.668, true, 3, 3, 3,
                                 true, true, true,
                                 types2, KIM_neigh_rvec_f, "box"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  Compute comp6(make_unique<Box>("B2", 2.668, false, 3, 3, 3,
                                 true, true, true,
                                 types2, KIM_neigh_rvec_f, "box"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  comp3.compute();
  comp4.compute();
  cout << comp3.get_energy_per_atom() << "  "
       << comp4.get_energy_per_atom() << endl;
  comp5.compute();
  comp6.compute();
  cout << comp5.get_energy_per_atom() << "  "
       << comp6.get_energy_per_atom() << endl;

  unique_ptr<Box> p = make_unique<Box>("B2", 2.668, false, 3, 3, 3,
                                       true, true, true,
                                       types2, KIM_neigh_rvec_f, "box");
  p->update_neighbor_list(3.0,0.0);
  p = comp1.change_box(move(p));
  comp1.compute();
  cout << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom" << endl;
  p = comp1.change_box(move(p));
  comp1.compute();
  cout << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom" << endl;

  const vector<int> typesSi{ 0 };
  const vector<int> typesC{ 1 };
  const vector<int> typesSiC{ 0, 1 };
  Compute comp(make_unique<Box>("B3", 4.359, false, 1, 1, 1,
                                true, true, true,
                                typesSiC, KIM_neigh_rvec_f, "box1"),
               "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  comp.compute();
  cout << comp.get_energy_per_atom() << " eV/atom" << endl;
  comp.move_atom(0, 0.1, 0.0, 0.0);
  vector<double> volumes1;
  vector<double> energies;
  const BMParams bmp = comp.bulk_modulus_energy(volumes1, energies, 0.05,
                                                false, false, true);
  cout << "E0  = " << bmp.E0 / comp.get_natoms() << " eV/atom\n"
       << "V0  = " << bmp.V0 / comp.get_natoms() << " A^3/atom\n"
       << "B0  = " << bmp.B0 * 160.2177 << " GPa\n"
       << "B0' = " << bmp.dB0_dp << endl;
  vector<double> volumes2;
  vector<double> pressures;
  const BMParams bmpp = comp.bulk_modulus_pressure(volumes2, pressures, 0.05,
                                                   false, false, true);
  cout << "V0  = " << bmpp.V0 / comp.get_natoms() << " A^3/atom\n"
       << "B0  = " << bmpp.B0 * 160.2177 << " GPa\n"
       << "B0' = " << bmpp.dB0_dp << endl;
  {
    ofstream outfile("EV.gnuplot");
    outfile << "set xzeroaxis" << endl
            << "plot "

            << bmp.E0<<" + 9*"<<bmpp.V0<<"*"<<bmpp.B0<<"/16 * "
            << "((("<<bmpp.V0<<"/x)**(2./3.) - 1)**3 * "<<bmpp.dB0_dp<<" + "
            << "(("<<bmpp.V0<<"/x)**(2./3.) - 1)**2 * "
            << "(6 - 4*("<<bmpp.V0<<"/x)**(2./3.)))"
            << "with lines lt 1 lc 1 title 'pressure-derived', "

            << bmp.E0<<" + 9*"<<bmp.V0<<"*"<<bmp.B0<<"/16 * "
            << "((("<<bmp.V0<<"/x)**(2./3.) - 1)**3 * "<<bmp.dB0_dp<<" + "
            << "(("<<bmp.V0<<"/x)**(2./3.) - 1)**2 * "
            << "(6 - 4*("<<bmp.V0<<"/x)**(2./3.)))"
            << "with lines lt 1 lc 3 title 'energy-derived', "

            << "'-' using 1:2 with points pt 7 title 'data'" << endl;
    for (unsigned i = 0; i != volumes1.size(); ++i)
      outfile << volumes1[i] <<  " " << energies[i] << endl;
    outfile << "e" << endl;
  }
  {
    ofstream outfile("pV.gnuplot");
    outfile << "set xzeroaxis" << endl
            << "plot "
            << "3*"<<bmpp.B0<<"/2 * "
            << "(("<<bmpp.V0<<"/x)**(7./3.) - ("<<bmpp.V0<<"/x)**(5./3.)) * "
            << "(1 + 0.75*("<<bmpp.dB0_dp<<" - 4) * "
            << "(("<<bmpp.V0<<"/x)**(2./3.) - 1))"
            << "with lines lt 1 lc 1 title 'pressure-derived', "
            << "3*"<<bmp.B0<<"/2 * "
            << "(("<<bmp.V0<<"/x)**(7./3.) - ("<<bmp.V0<<"/x)**(5./3.)) * "
            << "(1 + 0.75*("<<bmp.dB0_dp<<" - 4) * "
            << "(("<<bmp.V0<<"/x)**(2./3.) - 1))"
            << "with lines lt 1 lc 3 title 'energy-derived', "
            << "'-' using 1:2 with points pt 7 title 'data'" << endl;
    for (unsigned i = 0; i != volumes2.size(); ++i)
      outfile << volumes2[i] <<  " " << pressures[i] << endl;
    outfile << "e" << endl;
  }

  // Get elastic constants of silicon.
  Compute compSi(make_unique<Box>("diamond", 5.429, false, 1, 1, 1,
                                  true, true, true,
                                  typesSi, KIM_neigh_rvec_f, "boxxx"),
                 "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  comp.optimize_box(0.001, 10000);
  const double c11_0 = compSi.elastic_constant(1,1, 0.025, false);
  const double c11 = compSi.elastic_constant(1,1, 0.025, true);
  cout << "c11 = " << c11 * 160.2177 << "  " << c11_0 * 160.2177 << " GPa" << endl;
  const double c12_0 = compSi.elastic_constant(1,2, 0.025, false);
  const double c12 = compSi.elastic_constant(1,2, 0.025, true);
  cout << "c12 = " << c12 * 160.2177 << "  " << c12_0 * 160.2177 << " GPa" << endl;
  const double c44_0 = compSi.elastic_constant(4,4, 0.025, false);
  const double c44 = compSi.elastic_constant(4,4, 0.025, true);
  cout << "c44 = " << c44 * 160.2177 << "  " << c44_0 * 160.2177 << " GPa" << endl;

  // Let's check another model.
  const vector<int> typesAr = { 1 };
  Compute compAr(make_unique<Box>("fcc", 5.2, true, 1, 1, 1,
                                  true, true, true,
                                  typesAr, KIM_neigh_rvec_f, "box Ar"),
                 "ex_model_Ar_P_MLJ_NEIGH_RVEC_F");
  compAr.compute();
  cout << compAr.get_energy_per_atom() << endl;
  Voigt6<double> virial = compAr.get_virial();
  cout << virial.xx << "  "
       << virial.yy << "  "
       << virial.zz << "  "
       << virial.xy << "  "
       << virial.xz << "  "
       << virial.yz << endl;
  compAr.optimize_box(0.001, 100000, true);
  cout << "Lattice const "
       << compAr.get_a()[0] << "   "
       << compAr.get_b()[1] << "   "
       << compAr.get_c()[2] << endl;
  compAr.compute();
  cout << compAr.get_energy_per_atom() << endl;
  Voigt6<double> virial2 = compAr.get_virial();
  cout << virial2.xx << "  "
       << virial2.yy << "  "
       << virial2.zz << "  "
       << virial2.xy << "  "
       << virial2.xz << "  "
       << virial2.yz << endl;
  cout << "c11 = " << compAr.elastic_constant(1,1,0.005,true) * 160.2177 << endl;
  cout << "c22 = " << compAr.elastic_constant(2,2,0.005,true) * 160.2177 << endl;
  cout << "c33 = " << compAr.elastic_constant(3,3,0.005,true) * 160.2177 << endl;
  cout << "c44 = " << compAr.elastic_constant(4,4,0.005,true) * 160.2177 << endl;
  cout << "c55 = " << compAr.elastic_constant(5,5,0.005,true) * 160.2177 << endl;
  cout << "c66 = " << compAr.elastic_constant(6,6,0.005,true) * 160.2177 << endl;
  cout << "c12 = " << compAr.elastic_constant(1,2,0.005,true) * 160.2177 << endl;
  cout << "c13 = " << compAr.elastic_constant(1,3,0.005,true) * 160.2177 << endl;
  cout << "c14 = " << compAr.elastic_constant(1,4,0.005,true) * 160.2177 << endl;
  cout << "c15 = " << compAr.elastic_constant(1,5,0.005,true) * 160.2177 << endl;
  cout << "c16 = " << compAr.elastic_constant(1,6,0.005,true) * 160.2177 << endl;
  return 0;
  */
}
