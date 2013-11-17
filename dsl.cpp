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

#include "dsl.hpp"

#include <sstream>
#include <iterator>
#include <iostream>
#include <cctype>
#include <algorithm>

#include "utils.hpp"

using namespace std;
using namespace mytest;

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
  else if (to_lower(s) == "true")
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


void mytest::parse(const string& command,
                   map< string,unique_ptr<Box> >& boxes,
                   map<string,Compute>& computes) {
  vector<string> tokens = tokenize(command);
  if (tokens.size() == 0)
    // No command given.
    return;
  else if (tokens[0] == "box") {
    // Create a box
    if (tokens.size() < 13) {
      cout << "Not enough arguments." << endl;
      return;
    }
    try {
      vector<int> types;
      for (unsigned i = 12; i != tokens.size(); ++i)
        //types.psuh_back(compute.get_particle_type_code(tokens[i]));
        types.push_back(0); // TODO: cannot know the mapping here :-(    
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
  } else if (tokens[0] == "computer") {
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
  } else if (tokens[0] == "optimize_box")
    ;
  else if (tokens[0] == "optimize_positions")
    ;
  else if (tokens[0] == "bulk_modulus_energy")
    ;
  else if (tokens[0] == "bulk_modulus_pressure")
    ;
  else if (tokens[0] == "stiffness_tensor")
    ;
  else if (tokens[0] == "switch_boxes")
    ;
  else if (tokens[0] == "change_box")
    ;
  else
    cout << "Unknown command " << tokens[0] << endl;
}
