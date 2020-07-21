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

#include "utils.hpp"
#include "core.hpp"
#include "elastic.hpp"

using namespace std;
using namespace mytest;

void run_bench(const string& elem,
               const string& model) {
  random_device rd;
  mt19937 rng(rd());

  // Create simulation cells. //////////////////////////////////////////
  cout << "Creating box... " << flush;
  chrono::steady_clock::time_point start = chrono::steady_clock::now();
  vector<string> atomtypes = { elem };
  auto box = Box::random_box(50.0, 50.0, 50.0,
                             true, true, false,
                             2.1, atomtypes, "boxname", rng);
  chrono::steady_clock::time_point end = chrono::steady_clock::now();
  double duration = chrono::duration_cast<std::chrono::microseconds>(end - start).count(); // µs
  duration /= 1e6; // to s
  cout << "done in " << duration << "s." << endl;

  cout << "Repeating box... " << flush;
  start = chrono::steady_clock::now();
  box = box->repeat(3,1,1);
  end = chrono::steady_clock::now();
  duration = chrono::duration_cast<std::chrono::microseconds>(end - start).count(); // µs
  duration /= 1e6; // to s
  cout << "done in " << duration << "s." << endl;
  unsigned natoms = box->natoms;
  printf("natoms = %d\n", natoms);

  cout << "Creating neighbor lists... " << flush;
  start = chrono::steady_clock::now();
  Compute comp(make_unique<Box>(*box), model);
  end = chrono::steady_clock::now();
  duration = chrono::duration_cast<std::chrono::microseconds>(end - start).count(); // µs
  duration /= 1e6; // to s
  cout << "done in " << duration << "s." << endl;

  unsigned nsteps = 100;
  cout << "Computing " << nsteps << " times... " << flush;
  start = chrono::steady_clock::now();
  for (unsigned i = 0; i < nsteps; ++i) {
    comp.compute();
  }
  end = chrono::steady_clock::now();
  duration = chrono::duration_cast<std::chrono::microseconds>(end - start).count(); // µs
  duration /= 1e6; // to s
  cout << "done in " << duration << "s." << endl << endl;

  double d = duration / (natoms * nsteps);

  printf("This makes %g s/atom·step\n\n", d);
}


int main() {
  printf("Tersoff_LAMMPS_ErhartAlbe_2005_SiC__MO_903987585848_003\n");
  run_bench("Si", "Tersoff_LAMMPS_ErhartAlbe_2005_SiC__MO_903987585848_003");
  printf("\n");
  printf("Tersoff_LAMMPS_BjorkasNordlund_2007_Fe\n");
  run_bench("Fe", "Tersoff_LAMMPS_BjorkasNordlund_2007_Fe");
  return 0;
}
