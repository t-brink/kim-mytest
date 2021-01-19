/*
  Copyright (c) 2020 Tobias Brink

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
#include <chrono>
#include <sys/stat.h>

#include "utils.hpp"
#include "core.hpp"
#include "elastic.hpp"

using namespace std;
using namespace mytest;

int main() {
  random_device rd;
  mt19937 rng(rd());

  // Tersoff
  vector<string> atomtypes = { "Si" };
  auto box = Box::random_box(5.0, 5.0, 5.0,
                             true, true, false,
                             2.5, atomtypes, "boxname", rng);
  Compute comp(make_unique<Box>(*box), "Tersoff_LAMMPS_ErhartAlbe_2005_SiC__MO_903987585848_003");

  mkdir("param_output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  comp.write_kim_parameters("param_output/", "new_model_name");


  // Tersoff+ZBL
  vector<string> atomtypes_zbl = { "Fe" };
  auto box_zbl = Box::random_box(5.0, 5.0, 5.0,
                                 true, true, false,
                                 2.5, atomtypes_zbl, "boxname", rng);
  Compute comp_zbl(make_unique<Box>(*box_zbl), "Tersoff_LAMMPS_BjorkasNordlund_2007_Fe");

  mkdir("param_output_zbl", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  comp_zbl.write_kim_parameters("param_output_zbl", "new_model_name");


  return 0;
}
