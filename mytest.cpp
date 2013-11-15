#include <iostream>
#include <cstdio>

#include "mytestcore.hpp"
#include "utils.hpp"

using namespace std;
using namespace mytest;

int main() {
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
  p = comp1.change_box(move(p));
  comp1.compute();
  cout << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom" << endl;
  p = comp1.change_box(move(p));
  comp1.compute();
  cout << "comp1 = " << comp1.get_energy_per_atom() << " eV/atom" << endl;


  return 0;
}
