/*
  Copyright (c) 2013,2017 Tobias Brink

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

#include "dsl.hpp"

#include <sstream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <random>
#include <limits>

#include "utils.hpp"

using namespace std;
using namespace mytest;

namespace mytest {
std::random_device _DSL_RD_;
std::mt19937 _DSL_RNG_(_DSL_RD_());
}

vector<string> mytest::tokenize(const string& input) {
  istringstream iss(input);
  // Reading input from string stream separates at whitespace,
  // therefore this copies the whitespace-separated tokens to the
  // output vector.
  return { istream_iterator<string>{iss}, istream_iterator<string>{} };
}


static string to_lower(string data) {
  transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

static string to_upper(string data) {
  transform(data.begin(), data.end(), data.begin(), ::toupper);
  return data;
}

static double to_double(const string& s) {
  istringstream iss(s);
  double x;
  if (!(iss >> x))
    throw runtime_error("not a double: " + s);
  return x;
}

static unsigned to_unsigned(const string& s) {
  istringstream iss(s);
  unsigned x;
  if (!(iss >> x))
    throw runtime_error("not an unsigned int: " + s);
  return x;
}

static bool to_bool(const string& s) {
  if (to_lower(s) == "true")
    return true;
  else if (to_lower(s) == "false")
    return false;
  else
    throw runtime_error("not a boolean: " + s);
}

static KIMNeigh to_neighmode(const string& s) {
  if (to_upper(s) == "CLUSTER")
    return KIM_cluster;
  if (to_upper(s) == "MI_OPBC")
    return KIM_mi_opbc_f;
  if (to_upper(s) == "NEIGH_PURE")
    return KIM_neigh_pure_f;
  if (to_upper(s) == "NEIGH_RVEC")
    return KIM_neigh_rvec_f;
  else
    throw runtime_error("not a neighbor mode: " + s);
}


void mytest::parse(string command,
                   map< string,unique_ptr<Box> >& boxes,
                   map<string,Compute>& computes) {
  // Remove comments.
  auto comment_sign = find(command.begin(), command.end(), '#');
  command.erase(comment_sign, command.end());
  // Tokenize non-comment part.
  vector<string> tokens = tokenize(command);
  if (tokens.size() == 0) {
    // No command given.
    return;
  } else if (tokens[0] == "box") {
    // Create a box
    if (tokens.size() < 13) {
      cout << "Not enough arguments." << endl;
      return;
    }
    try {
      vector<string> types;
      for (unsigned i = 12; i != tokens.size(); ++i)
        types.push_back(tokens[i]);
      auto p = make_unique<Box>(tokens[2], to_double(tokens[3]), // Lattice
                                to_bool(tokens[4]), // Cubic?
                                to_unsigned(tokens[5]), // repeat a
                                to_unsigned(tokens[6]), // repeat b
                                to_unsigned(tokens[7]), // repeat c
                                to_bool(tokens[8]),  // periodic a
                                to_bool(tokens[9]),  // periodic b
                                to_bool(tokens[10]), // periodic c
                                types, // Atom types
                                to_neighmode(tokens[11]), // Neighbor list mode
                                tokens[1] // Name
                                );
      boxes[tokens[1]] = move(p);
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "random_box") {
    // Create a box
    if (tokens.size() < 8) {
      cout << "Not enough arguments." << endl;
      return;
    }
    try {
      unique_ptr<Box> p =
        Box::random_box(10.0, 10.0, 10.0,
                        to_bool(tokens[2]), // periodic a
                        to_bool(tokens[3]), // periodic b
                        to_bool(tokens[4]), // periodic c
                        to_double(tokens[5]), // min_dist between atoms
                        tokens[7], // Atom type (TODO: more than one)    
                        to_neighmode(tokens[6]), // Neighbor list mode
                        tokens[1], // name
                        _DSL_RNG_);
      boxes[tokens[1]] = move(p);
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "copy_box") {
    // Create a box
    if (tokens.size() != 3) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    auto box = boxes.find(tokens[1]);
    if (box == boxes.end()) {
      cout << "Unknown box: " << tokens[1] << endl;
      return;
    }
    try {
      auto p = make_unique<Box>(*(box->second), tokens[2]);
      boxes[tokens[2]] = move(p);
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "delete_atom") {
    if (tokens.size() < 2) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 3) {
      cout << "Too many arguments." << endl;
      return;
    }
    auto box = boxes.find(tokens[1]);
    if (box == boxes.end()) {
      cout << "Unknown box: " << tokens[1] << endl;
      return;
    }
    unsigned i;
    if (tokens.size() == 2) {
      // Delete random atom
      uniform_int_distribution<unsigned> uni_rng(0, box->second->natoms - 1);
      i = uni_rng(_DSL_RNG_);
    } else {
      i = to_unsigned(tokens[2]);
      if (i >= box->second->natoms) {
        cout << "Atom " << i << " does not exist." << endl;
        return;
      }
    }
    boxes[box->first] = box->second->delete_atom(i);
  } else if (tokens[0] == "model") {
    // Create a compute
    if (tokens.size() < 4) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 4) {
      cout << "Too many arguments." << endl;
      return;
    }
    try {
      const auto iter = boxes.find(tokens[2]);
      if (iter == boxes.end()) {
        cout << "Unknown box: " << tokens[2] << endl;
        return;
      }
      auto b = make_unique<Box>(*(iter->second));
      Compute c(move(b), tokens[3]);
      computes.emplace(piecewise_construct,
                       forward_as_tuple(tokens[1]),
                       forward_as_tuple(make_unique<Box>(*(iter->second)),
                                        tokens[3]));
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "deform_box") {
    if (tokens.size() < 8) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 8) {
      cout << "Too many arguments." << endl;
      return;
    }
    auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown compute: " << tokens[1] << endl;
      return;
    }
    iter->second.deform(Voigt6<double>(to_double(tokens[2]),
                                       to_double(tokens[3]),
                                       to_double(tokens[4]),
                                       to_double(tokens[5]),
                                       to_double(tokens[6]),
                                       to_double(tokens[7])));
  } else if (tokens[0] == "write_box") {
    if (tokens.size() < 3) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 3) {
      cout << "Too many arguments." << endl;
      return;
    }
    auto iter = boxes.find(tokens[1]);
    if (iter == boxes.end()) {
      cout << "Unknown box: " << tokens[1] << endl;
      return;
    }
    iter->second->write_to(tokens[2]);
  } else if (tokens[0] == "compute") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    try {
      iter->second.compute();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    cout << "Cohesive energy: "
         << iter->second.get_energy_per_atom() << " eV/atom" << endl;
  } else if (tokens[0] == "check_testfile") {
    if (tokens.size() != 4) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    // Read the box data.
    ifstream testfile(tokens[2]); // TODO: error checking for I/O   
    string comment;
    getline(testfile, comment); // first line is comment line.
    unsigned natoms;
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    bool pa, pb, pc;
    testfile >> natoms
             >> ax >> ay >> az >> pa
             >> bx >> by >> bz >> pb
             >> cx >> cy >> cz >> pc;
    auto coordinates = make_unique< Array2D<double> >(natoms, 3);
    auto types = make_unique< Array1DInit<string> >(natoms);
    Array1D<double> energies(natoms);
    Array2D<double> forces(natoms, 3);
    Array2D<double> stresses(natoms, 6);
    for (unsigned i = 0; i < natoms; ++i) {
      testfile >> (*types)(i)
               >> (*coordinates)(i, 0)
               >> (*coordinates)(i, 1)
               >> (*coordinates)(i, 2)
               >> energies(i)
               >> forces(i, 0)
               >> forces(i, 1)
               >> forces(i, 2)
               >> stresses(i, 0)
               >> stresses(i, 1)
               >> stresses(i, 2)
               >> stresses(i, 3)
               >> stresses(i, 4)
               >> stresses(i, 5);
    }
    // Make a box.
    auto box = make_unique<Box>(Vec3D<double>(ax, ay, az),
                                Vec3D<double>(bx, by, bz),
                                Vec3D<double>(cx, cy, cz),
                                pa, pb, pc,
                                move(coordinates),
                                move(types),
                                KIM_neigh_pure_f,
                                tokens[3]);
    // Assign box to compute.
    iter->second.change_box(move(make_unique<Box>(*box)));
    // Store box under its name.
    boxes[tokens[3]] = move(box);
    // Compute.
    try {
      iter->second.compute();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    double max_err_energy = 0.0;
    double max_err_force = 0.0;
    double max_err_stress = 0.0;
    // Compare.
    for (unsigned i = 0; i < natoms; ++i) {
      //
      double delta_E = abs(iter->second.get_energy(i) - energies(i));
      if (delta_E > max_err_energy) max_err_energy = delta_E;
      //
      for (unsigned dim = 0; dim < 3; ++dim) {
        double delta_F = abs(iter->second.get_force(i,dim) - forces(i,dim));
        if (delta_F > max_err_force) max_err_force = delta_F;
      }
      //
      Voigt6<double> v = iter->second.get_virial(i);
      for (unsigned dim = 0; dim < 6; ++dim) {
        double delta_v = abs(v(dim) - stresses(i,dim));
        if (delta_v > max_err_stress) max_err_stress = delta_v;
      }
    }
    // Output.
    cout << comment << endl;
    bool has_err = (abs(max_err_energy) > 1e-5)
                || (abs(max_err_force) > 1e-5);
    bool has_err_stress = abs(max_err_stress) > 1e-4;
    cout << "Maximum particle energy deviation: " << max_err_energy;
    if (has_err)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    else if (has_err_stress)
      cout << "     !!!"; // just for grepping
    cout << endl;
    //
    cout << "Maximum force deviation: " << max_err_force;
    if (has_err)
      cout << "               !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    else if (has_err_stress)
      cout << "     !!!"; // just for grepping
    cout << endl;
    //
    cout << "Maximum particle virial deviation: " << max_err_stress;
    if (has_err_stress)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
    //
  } else if (tokens[0] == "optimize_box") {
    // Optimize box.
    if (tokens.size() < 2) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 3) {
      cout << "Too many arguments." << endl;
      return;
    }
    // Do it.
    try {
      auto iter = computes.find(tokens[1]);
      if (iter == computes.end()) {
        cout << "Unknown model: " << tokens[1] << endl;
        return;
      }
      Compute& comp = iter->second;
      if (tokens.size() == 2)
        comp.optimize_box(0.001, 10000, false);
      else
        comp.optimize_box(0.001, 10000, to_bool(tokens[2]));
      cout << "New box vectors after " << comp.get_optimization_steps() << " steps:"
           << endl;
      printf("a: %10.4g %10.4g %10.4g\n",comp.get_a()[0],comp.get_a()[1],comp.get_a()[2]);
      printf("b: %10.4g %10.4g %10.4g\n",comp.get_b()[0],comp.get_b()[1],comp.get_b()[2]);
      printf("c: %10.4g %10.4g %10.4g\n",comp.get_c()[0],comp.get_c()[1],comp.get_c()[2]);
      cout << "Cohesive energy is now " << comp.get_energy_per_atom() << " eV/atom"
           << endl;
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "optimize_positions") {
    // Optimize positions.
    if (tokens.size() < 2) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 2) {
      cout << "Too many arguments." << endl;
      return;
    }
    // Do it.
    try {
      auto iter = computes.find(tokens[1]);
      if (iter == computes.end()) {
        cout << "Unknown model: " << tokens[1] << endl;
        return;
      }
      Compute& comp = iter->second;
      comp.optimize_positions(0.001, 10000);
      cout << "After " << comp.get_optimization_steps() << " steps:"
           << endl;
      cout << "Cohesive energy is now " << comp.get_energy_per_atom() << " eV/atom"
           << endl;
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "bulk_modulus_energy") {
    if (tokens.size() < 3) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 9) {
      cout << "Too many arguments." << endl;
      return;
    }
    const auto comp = computes.find(tokens[1]);
    if (comp == computes.end()) {
      cout << "Unknown model: " << tokens[1] << endl;
      return;
    }
    const bool c_to_a = (tokens.size() > 3) ? to_bool(tokens[3])
                                             : false;
    const bool b_to_a = (tokens.size() > 4) ? to_bool(tokens[4])
                                             : false;
    const bool positions = (tokens.size() > 5) ? to_bool(tokens[5])
                                                : false;
    const bool angle_ab = (tokens.size() > 6) ? to_bool(tokens[6])
                                               : false;
    const bool angle_ac = (tokens.size() > 7) ? to_bool(tokens[7])
                                               : false;
    const bool angle_bc = (tokens.size() > 8) ? to_bool(tokens[8])
                                               : false;
    const BMParams bmp = comp->second.bulk_modulus_energy(to_double(tokens[2]),
                                                          c_to_a, b_to_a,
                                                          positions,
                                                          angle_ab, angle_ac,
                                                          angle_bc);
    const unsigned natoms = comp->second.get_natoms();
    cout << "E0  = " << bmp.E0 / natoms << " eV/atom" << endl;
    cout << "V0  = " << bmp.V0 / natoms << " Ang/atom" << endl;
    cout << "B0  = " << bmp.B0 * 160.2177 << " GPa" << endl;
    cout << "B0' = " << bmp.dB0_dp << endl;

  } else if (tokens[0] == "bulk_modulus_pressure") {
    if (tokens.size() < 3) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 9) {
      cout << "Too many arguments." << endl;
      return;
    }
    const auto comp = computes.find(tokens[1]);
    if (comp == computes.end()) {
      cout << "Unknown model: " << tokens[1] << endl;
      return;
    }
    const bool c_to_a = (tokens.size() > 3) ? to_bool(tokens[3])
                                             : false;
    const bool b_to_a = (tokens.size() > 4) ? to_bool(tokens[4])
                                             : false;
    const bool positions = (tokens.size() > 5) ? to_bool(tokens[5])
                                                : false;
    const bool angle_ab = (tokens.size() > 6) ? to_bool(tokens[6])
                                               : false;
    const bool angle_ac = (tokens.size() > 7) ? to_bool(tokens[7])
                                               : false;
    const bool angle_bc = (tokens.size() > 8) ? to_bool(tokens[8])
                                               : false;
    const BMParams bmp =comp->second.bulk_modulus_pressure(to_double(tokens[2]),
                                                           c_to_a, b_to_a,
                                                           positions,
                                                           angle_ab, angle_ac,
                                                           angle_bc);
    const unsigned natoms = comp->second.get_natoms();
    cout << "V0  = " << bmp.V0 / natoms << " Ang/atom" << endl;
    cout << "B0  = " << bmp.B0 * 160.2177 << " GPa" << endl;
    cout << "B0' = " << bmp.dB0_dp << endl;
  } else if (tokens[0] == "stiffness_tensor") {
    if (tokens.size() < 3) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 4) {
      cout << "Too many arguments." << endl;
      return;
    }
    const auto comp = computes.find(tokens[1]);
    if (comp == computes.end()) {
      cout << "Unknown model: " << tokens[1] << endl;
      return;
    }
    const double max_strain = to_double(tokens[2]);
    const bool positions = (tokens.size() > 3) ? to_bool(tokens[3])
                                                : false;
    for (unsigned i = 1; i <= 6; ++i) {
      for (unsigned j = 1; j <= 6; ++j) {
        const double cij = comp->second.elastic_constant(i,j,max_strain,
                                                         positions);
        printf("%8.2f  ", cij * 160.2177);
      }
      cout << "GPa" << endl;
    }
  } else if (tokens[0] == "switch_boxes")
    // TODO?                                                 
    ;
  else if (tokens[0] == "change_box") {
    // Change box in Compute.
    if (tokens.size() < 3) {
      cout << "Not enough arguments." << endl;
      return;
    } else if (tokens.size() > 3) {
      cout << "Too many arguments." << endl;
      return;
    }
    // Do it.
    try {
      const auto iter = boxes.find(tokens[2]);
      if (iter == boxes.end()) {
        cout << "Unknown box: " << tokens[2] << endl;
        return;
      }
      auto b = make_unique<Box>(*(iter->second));
      const auto iter_comp = computes.find(tokens[1]);
      if (iter_comp == computes.end()) {
        cout << "Unknown model: " << tokens[1] << endl;
        return;
      }
      iter_comp->second.change_box(move(b));
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
  } else if (tokens[0] == "scale") {
    if (tokens.size() != 3 && tokens.size() != 5) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto comp = computes.find(tokens[1]);
    if (comp == computes.end()) {
      cout << "Unknown model: " << tokens[1] << endl;
      return;
    }
    if (tokens.size() == 3) {
      cout << "Not implemented..." << endl;
      return;
    } else {
      comp->second.scale_box(to_double(tokens[2]),
                             to_double(tokens[3]),
                             to_double(tokens[4]));
    }
  } else if (tokens[0] == "numer_forces_deriv") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    double max_diff_force = 666;
    try {
      max_diff_force = iter->second.max_diff_force();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    cout << "Maximum force deviation: " << max_diff_force;
    if (abs(max_diff_force) > 1e-5)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
    /*
    cout << endl;
    iter->second.compute(); // otherwise we get the displaced one!!
    for (unsigned i = 0; i < iter->second.get_natoms(); ++i) {
      printf("Forces on atom %3d:", i);
      for (unsigned dim = 0; dim < 3; ++dim)
        //printf(" %+10.3e", iter->second.get_force(i, dim));
        printf(" %+16.10f", iter->second.get_force(i, dim));
      printf("   pos:");
      Vec3D<double> pos = iter->second.get_position(i);
      for (unsigned dim = 0; dim < 3; ++dim)
        printf(" %7.4f", pos[dim]);
      printf("     %7.4f", iter->second.get_energy(i));
      printf("\n");
    }
    */
  } else if (tokens[0] == "numer_forces_deriv_fast") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    double max_diff_force = 666;
    try {
      max_diff_force = iter->second.max_diff_force_fast();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    cout << "Maximum force deviation: " << max_diff_force;
    if (abs(max_diff_force) > 1e-5)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
  } else if (tokens[0] == "diff_total_energy_vs_particle_energy") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    try {
      iter->second.compute();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    unsigned N = iter->second.get_natoms();
    double E = iter->second.get_energy();
    double Epart = 0.0;
    for (unsigned i = 0; i < N; ++i) {
      Epart += iter->second.get_energy(i);
    }
    double diff = abs(E - Epart);
    cout << "Energy diff: " << diff;
    if (diff > 1e-5)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
  } else if (tokens[0] == "diff_total_virial_vs_particle_virial") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    try {
      iter->second.compute();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    unsigned N = iter->second.get_natoms();
    Voigt6<double> virial = iter->second.get_virial();
    double vxx = 0.0;
    double vyy = 0.0;
    double vzz = 0.0;
    double vyz = 0.0;
    double vxz = 0.0;
    double vxy = 0.0;
    for (unsigned i = 0; i < N; ++i) {
      Voigt6<double> vi = iter->second.get_virial(i);
      vxx += vi.xx;
      vyy += vi.yy;
      vzz += vi.zz;
      vyz += vi.yz;
      vxz += vi.xz;
      vxy += vi.xy;
    }
    double
      dxx = vxx - virial.xx,
      dxy = vxy - virial.xy,
      dxz = vxz - virial.xz,
      dyy = vyy - virial.yy,
      dyz = vyz - virial.yz,
      dzz = vzz - virial.zz;
    bool too_large = (abs(dxx) > 1e-5)
                     || (abs(dxy) > 1e-5)
                     || (abs(dxz) > 1e-5)
                     || (abs(dyy) > 1e-5)
                     || (abs(dyz) > 1e-5)
                     || (abs(dzz) > 1e-5);
    cout << "Virial diff:";
    if (too_large)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
    printf("  %+16.10f %+16.10f %+16.10f", dxx, dxy, dxz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
    printf("  %16s %+16.10f %+16.10f", "", dyy, dyz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
    printf("  %16s %16s %+16.10f", "", "", dzz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
  } else if (tokens[0] == "diff_total_virial_vs_virial_from_forces") {
    if (tokens.size() != 2) {
      cout << "Wrong number of arguments." << endl;
      return;
    }
    const auto iter = computes.find(tokens[1]);
    if (iter == computes.end()) {
      cout << "Unknown computer: " << tokens[1] << endl;
      return;
    }
    try {
      iter->second.compute();
    } catch (const exception& e) {
      cout << e.what() << endl;
      return;
    }
    Voigt6<double> virial = iter->second.get_virial();
    Voigt6<double> virial_forces = iter->second.get_virial_from_forces();
    double
      dxx = virial_forces.xx - virial.xx,
      dxy = virial_forces.xy - virial.xy,
      dxz = virial_forces.xz - virial.xz,
      dyy = virial_forces.yy - virial.yy,
      dyz = virial_forces.yz - virial.yz,
      dzz = virial_forces.zz - virial.zz;
    bool too_large = (abs(dxx) > 1e-5)
                     || (abs(dxy) > 1e-5)
                     || (abs(dxz) > 1e-5)
                     || (abs(dyy) > 1e-5)
                     || (abs(dyz) > 1e-5)
                     || (abs(dzz) > 1e-5);
    cout << "Virial diff:";
    if (too_large)
      cout << "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    cout << endl;
    printf("  %+16.10f %+16.10f %+16.10f", dxx, dxy, dxz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
    printf("  %16s %+16.10f %+16.10f", "", dyy, dyz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
    printf("  %16s %16s %+16.10f", "", "", dzz);
    if (too_large) printf("     !!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    printf("\n");
    /*
    cout << "Virial model:" << endl;
    printf("  %+16.10f %+16.10f %+16.10f\n",
           virial.xx,
           virial.xy,
           virial.xz);
    printf("  %16s %+16.10f %+16.10f\n",
           "",
           virial.yy,
           virial.yz);
    printf("  %16s %16s %+16.10f\n",
           "",
           "",
           virial.zz);
    cout << "Virial forces:" << endl;
    printf("  %+16.10f %+16.10f %+16.10f\n",
           virial_forces.xx,
           virial_forces.xy,
           virial_forces.xz);
    printf("  %16s %+16.10f %+16.10f\n",
           "",
           virial_forces.yy,
           virial_forces.yz);
    printf("  %16s %16s %+16.10f\n",
           "",
           "",
           virial_forces.zz);
    */
  } else if (tokens[0] == "print") {
    if (command.size() > 6)
      cout << command.substr(6) << flush;
  } else if (tokens[0] == "println") {
    if (command.size() > 8)
      cout << command.substr(8);
    cout << endl;
  } else
    cout << "Unknown command " << tokens[0] << endl;
}
