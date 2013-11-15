#include <iostream>
#include <cstdio>

#include "mytestcore.hpp"
#include "utils.hpp"

using namespace std;
using namespace mytest;

int main() {
  const vector<int> types{ 0 };
  Compute comp1(make_unique<Box>("fcc", 3.940, true, 3, 3, 3,
                                 true, true, true,
                                 types, KIM_neigh_rvec_f, "box1"),
                "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
  Compute comp2(make_unique<Box>("sc", 2.525, true, 3, 3, 3,
                                 true, true, true,
                                 types, KIM_neigh_rvec_f, "box2"),
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

  return 0;
}
