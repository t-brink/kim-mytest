/*
  Copyright (c) 2013,2014,2017,2018,2019 Tobias Brink

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

//#include <iostream>
//#include <cstdio>
#include <utility>
#include <cmath>
#include <iostream>
//#include <chrono>

#include <nlopt.hpp>

#include "core.hpp"

using namespace std;
using namespace mytest;

const double DELTA = 1e-10; // For double comparison.

// Box /////////////////////////////////////////////////////////////////

Box::Box(const Vec3D<double>& a_in,
         const Vec3D<double>& b_in,
         const Vec3D<double>& c_in,
         bool periodic_a, bool periodic_b, bool periodic_c,
         unique_ptr< Array2D<double> > coordinates,
         unique_ptr< Array1DInit<string> > types_in,
         const string& name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(move(*coordinates)), // Public data.
    types(move(*types_in)),
    box_side_lengths_(a_in.abs(), b_in.abs(), c_in.abs()), // Private data.
    a_(a_in), b_(b_in), c_(c_in),
    periodic_(periodic_a, periodic_b, periodic_c),
    natoms_(coordinates->extent(0)), nghosts_(0), nall_(natoms_),
    name_(name),
    neigh_list_(natoms_),
    contributing_(natoms_, 1) // init to ones for all "real" atoms
{
  if (positions.extent(0) != types.extent(0))
    throw runtime_error("types must have the same length as the "
                        "first dimension of positions.");
  if (positions.extent(1) != 3)
    throw runtime_error("positions must be a n*3 array.");
}

Box::Box(const std::string& lattice, double lattice_const, bool cubic,
         unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
         bool periodic_a, bool periodic_b, bool periodic_c,
         const vector<string>& types_in,
         const std::string& name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(atoms_per_unit_cell(lattice, cubic) // Public data.
              * repeat_a * repeat_b * repeat_c, 3),
    types(positions.extent(0)),
    periodic_(periodic_a, periodic_b, periodic_c), // Private data.
    natoms_(positions.extent(0)), nghosts_(0), nall_(natoms_),
    name_(name),
    contributing_(natoms_, 1) // init to ones for all "real" atoms
{
  if (types_in.size() != species_per_unit_cell(lattice))
    throw runtime_error("Wrong number of types for this lattice.");

  // The unit cell vectors.
  Vec3D<double> unit_a;
  Vec3D<double> unit_b;
  Vec3D<double> unit_c;
  typedef pair<Vec3D<double>,string> t_atom;
  vector<t_atom> atoms;   // (position, type)

  if (lattice == "dimer") {
    // Arbitrary big box size.
    unit_a = Vec3D<double>(100*lattice_const, 0.0, 0.0);
    unit_b = Vec3D<double>(0.0, 100*lattice_const, 0.0);
    unit_c = Vec3D<double>(0.0, 0.0, 100*lattice_const);
    atoms.push_back(t_atom(Vec3D<double>(49.5*lattice_const,
                                         50.0*lattice_const,
                                         50.0*lattice_const),
                           types_in[0]));
    atoms.push_back(t_atom(Vec3D<double>(50.5*lattice_const,
                                         50.0*lattice_const,
                                         50.0*lattice_const),
                           types_in[1]));
  } else if (lattice == "sc") {
    unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
    unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
    unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
    atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0), types_in[0]));
  } else if (lattice == "bcc") {
    if (cubic) {
      unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
      unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
      unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0),
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.5, 0.5) * lattice_const,
                             types_in[0]));
    } else {
      unit_a = Vec3D<double>(-0.5, 0.5, 0.5) * lattice_const;
      unit_b = Vec3D<double>(0.5, -0.5, 0.5) * lattice_const;
      unit_c = Vec3D<double>(0.5, 0.5, -0.5) * lattice_const;
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0), types_in[0]));
    }
  } else if (lattice == "fcc") {
    if (cubic) {
      unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
      unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
      unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0),
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.5, 0.5) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.0, 0.5) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.5, 0.0) * lattice_const,
                             types_in[0]));
    } else {
      unit_a = Vec3D<double>(0.0, 0.5, 0.5) * lattice_const;
      unit_b = Vec3D<double>(0.5, 0.0, 0.5) * lattice_const;
      unit_c = Vec3D<double>(0.5, 0.5, 0.0) * lattice_const;
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0), types_in[0]));
    }
  } else if (lattice == "diamond") {
    if (cubic) {
      unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
      unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
      unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.00, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.25, 0.25) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.50, 0.50) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.75, 0.75) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.50, 0.00, 0.50) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.75, 0.25, 0.75) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.50, 0.50, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.75, 0.75, 0.25) * lattice_const,
                             types_in[0]));

    } else {
      unit_a = Vec3D<double>(0.0, 0.5, 0.5) * lattice_const;
      unit_b = Vec3D<double>(0.5, 0.0, 0.5) * lattice_const;
      unit_c = Vec3D<double>(0.5, 0.5, 0.0) * lattice_const;
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.00, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.25, 0.25) * lattice_const,
                             types_in[0]));
    }
  } else if (lattice == "B1") { // NaCl
    if (cubic) {
      unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
      unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
      unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.0, 0.0) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.5, 0.5) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.5, 0.5) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.0, 0.5) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.5) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.5, 0.0) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.5, 0.0) * lattice_const,
                             types_in[1]));

    } else {
      unit_a = Vec3D<double>(0.0, 0.5, 0.5) * lattice_const;
      unit_b = Vec3D<double>(0.5, 0.0, 0.5) * lattice_const;
      unit_c = Vec3D<double>(0.5, 0.5, 0.0) * lattice_const;
      atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.5, 0.0, 0.0) * lattice_const,
                             types_in[1]));
    }
  } else if (lattice == "B2") { // CsCl
    unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
    unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
    unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
    atoms.push_back(t_atom(Vec3D<double>(0.0, 0.0, 0.0) * lattice_const,
                           types_in[0]));
    atoms.push_back(t_atom(Vec3D<double>(0.5, 0.5, 0.5) * lattice_const,
                           types_in[1]));
  } else if (lattice == "B3") { // zincblende
    if (cubic) {
      unit_a = Vec3D<double>(lattice_const, 0.0, 0.0);
      unit_b = Vec3D<double>(0.0, lattice_const, 0.0);
      unit_c = Vec3D<double>(0.0, 0.0, lattice_const);
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.00, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.25, 0.25) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.50, 0.50) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.75, 0.75) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.50, 0.00, 0.50) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.75, 0.25, 0.75) * lattice_const,
                             types_in[1]));
      atoms.push_back(t_atom(Vec3D<double>(0.50, 0.50, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.75, 0.75, 0.25) * lattice_const,
                             types_in[1]));

    } else {
      unit_a = Vec3D<double>(0.0, 0.5, 0.5) * lattice_const;
      unit_b = Vec3D<double>(0.5, 0.0, 0.5) * lattice_const;
      unit_c = Vec3D<double>(0.5, 0.5, 0.0) * lattice_const;
      atoms.push_back(t_atom(Vec3D<double>(0.00, 0.00, 0.00) * lattice_const,
                             types_in[0]));
      atoms.push_back(t_atom(Vec3D<double>(0.25, 0.25, 0.25) * lattice_const,
                             types_in[1]));
    }
  } else {
    throw runtime_error("unknown lattice");
  }
  // Create box.
  a_ = double(repeat_a) * unit_a; box_side_lengths_[0] = a_.abs();
  b_ = double(repeat_b) * unit_b; box_side_lengths_[1] = b_.abs();
  c_ = double(repeat_c) * unit_c; box_side_lengths_[2] = c_.abs();
  // Fill atoms.
  unsigned i = 0;
  for (unsigned aa = 0; aa != repeat_a; ++aa)
    for (unsigned bb = 0; bb != repeat_b; ++bb)
      for (unsigned cc = 0; cc != repeat_c; ++cc) {
        // Get offset.
        const Vec3D<double> offset =
          double(aa)*unit_a + double(bb)*unit_b + double(cc)*unit_c;
        // Assign atoms.
        for (const t_atom& atom : atoms) {
          for (unsigned j = 0; j != 3; ++j)
            positions(i, j) = atom.first[j] + offset[j];
          types(i) = atom.second;
          ++i;
        }
      }
  if (i != natoms_)
    throw runtime_error("BUG in the implementation, not all atoms are filled.");

  // Reserve neighbor list space.
  neigh_list_.resize(natoms_);
}

Box::Box(const Box& other, const std::string& new_name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(other.positions.extent(0), // Public data.
              other.positions.extent(1)),
    types(other.types.extent(0)),
    box_side_lengths_(other.box_side_lengths_), // Private data.
    a_(other.a_), b_(other.b_), c_(other.c_),
    periodic_(other.periodic_),
    natoms_(other.natoms_), nghosts_(other.nghosts_), nall_(other.nall_),
    name_(new_name),
    ghost_shells(other.ghost_shells),
    /*
    ghost_positions(make_unique<Array2D<double>>(
                                    other.ghost_positions->extent(0),
                                    other.ghost_positions->extent(1))),
    ghost_types(make_unique<Array1D<int>>(other.ghost_types->extent(0))),
    */
    neigh_list_(other.neigh_list_),
    contributing_(other.contributing_)
{
  for (int i = 0; i != positions.extent(0); ++i)
    for (int j = 0; j != positions.extent(1); ++j)
      positions(i,j) = other.positions(i,j);
  for (int i = 0; i != types.extent(0); ++i)
    types(i) = other.types(i);
  if (other.ghost_positions && other.ghost_types) {
    ghost_positions = make_unique<Array2D<double>>(other.ghost_positions->extent(0),
                                                   other.ghost_positions->extent(1));
    ghost_types = make_unique<Array1D<int>>(other.ghost_types->extent(0));
    for (int i = 0; i != ghost_positions->extent(0); ++i)
      for (int j = 0; j != ghost_positions->extent(1); ++j)
        (*ghost_positions)(i,j) = (*other.ghost_positions)(i,j);
    for (int i = 0; i != ghost_types->extent(0); ++i)
      (*ghost_types)(i) = (*other.ghost_types)(i);
  }
}

Box::Box(const Box& other)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(other.positions.extent(0), // Public data.
              other.positions.extent(1)),
    types(other.types.extent(0)),
    box_side_lengths_(other.box_side_lengths_), // Private data.
    a_(other.a_), b_(other.b_), c_(other.c_),
    periodic_(other.periodic_),
    natoms_(other.natoms_), nghosts_(other.nghosts_), nall_(other.nall_),
    name_(other.name_),
    ghost_shells(other.ghost_shells),
    /*
    ghost_positions(make_unique<Array2D<double>>(
                                    other.ghost_positions->extent(0),
                                    other.ghost_positions->extent(1))),
    ghost_types(make_unique<Array1D<int>>(other.ghost_types->extent(0))),
    */
    neigh_list_(other.neigh_list_),
    contributing_(other.contributing_)
{
  for (int i = 0; i != positions.extent(0); ++i)
    for (int j = 0; j != positions.extent(1); ++j)
      positions(i,j) = other.positions(i,j);
  for (int i = 0; i != types.extent(0); ++i)
    types(i) = other.types(i);
  if (other.ghost_positions && other.ghost_types) {
    ghost_positions = make_unique<Array2D<double>>(other.ghost_positions->extent(0),
                                                   other.ghost_positions->extent(1));
    ghost_types = make_unique<Array1D<int>>(other.ghost_types->extent(0));
    for (int i = 0; i != ghost_positions->extent(0); ++i)
      for (int j = 0; j != ghost_positions->extent(1); ++j)
        (*ghost_positions)(i,j) = (*other.ghost_positions)(i,j);
    for (int i = 0; i != ghost_types->extent(0); ++i)
      (*ghost_types)(i) = (*other.ghost_types)(i);
  }
}



bool Box::update_neighbor_list(double cutoff, double skin,
                               const map<string,int>& typemap) {
  const double cut = (1 + skin) * cutoff;
  ghost_shells = calc_number_of_ghost_shells(cut);
  unsigned new_nghosts;
  // We store any ghost atoms that we may need.
  new_nghosts = natoms_ * ((2*ghost_shells[0]+1)
                           * (2*ghost_shells[1]+1)
                           * (2*ghost_shells[2]+1) - 1);
  // If needed, resize ghost arrays. If the arrays are not yet
  // allocated, allocate them (ghost_position == nullptr).
  bool reallocated = false;
  if (new_nghosts != nghosts_ || !ghost_positions) {
    // Number of ghosts changed, re-allocate memory.
    nghosts_ = new_nghosts;
    nall_ = natoms_ + nghosts_;
    ghost_positions = make_unique< Array2D<double> >(nall_, 3);
    ghost_types = make_unique< Array1D<int> >(nall_);
    contributing_.resize(nall_, 0); // Everything that the vector is
                                    // resized to consists only of
                                    // ghost atoms, thus init to zero.
    reallocated = true;
  }

  // Update ghosts.
  update_ghosts(typemap);

  // Clear neighbor lists.
  for (auto& l : neigh_list_)
    l.clear();

  // Fill neighbor lists again.
  const double cutsq = cut * cut;
  for (unsigned i = 0; i != natoms_; ++i) {
    // Go over neighbors in the central cell.
    for (unsigned j = i+1; j != natoms_; ++j) {
      const double dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
      const double dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
      const double dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
      if (dx*dx + dy*dy + dz*dz < cutsq) {
        neigh_list_[i].push_back(j);
        neigh_list_[j].push_back(i);
      }
    }
    // Iterate over ghosts.
    for (unsigned j = natoms_; j != nall_; ++j) {
      const double dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
      const double dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
      const double dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
      if (dx*dx + dy*dy + dz*dz < cutsq) {
        neigh_list_[i].push_back(j);
      }
    }
  }

  return reallocated;
}

void Box::update_ghosts(const map<string,int>& typemap) {
  // Copy original atoms, do not wrap for PBCs.
  for (unsigned i = 0; i != natoms_; ++i) {
    (*ghost_positions)(i, 0) = positions(i, 0);
    (*ghost_positions)(i, 1) = positions(i, 1);
    (*ghost_positions)(i, 2) = positions(i, 2);
    (*ghost_types)(i) = typemap.at(types(i));
  }
  unsigned ii = natoms_;
  const int alo = -static_cast<int>(ghost_shells[0]);
  const int ahi = ghost_shells[0];
  const int blo = -static_cast<int>(ghost_shells[1]);
  const int bhi = ghost_shells[1];
  const int clo = -static_cast<int>(ghost_shells[2]);
  const int chi = ghost_shells[2];
  for (int aa = alo; aa <= ahi; ++aa)
    for (int bb = blo; bb <= bhi; ++bb)
      for (int cc = clo; cc <= chi; ++cc) {
        if (aa == 0 && bb == 0 && cc == 0) continue;
        const Vec3D<double> offset =
          double(aa)*a + double(bb)*b + double(cc)*c;
        for (unsigned i = 0; i != natoms_; ++i) {
          (*ghost_positions)(ii, 0) = (*ghost_positions)(i, 0) + offset[0];
          (*ghost_positions)(ii, 1) = (*ghost_positions)(i, 1) + offset[1];
          (*ghost_positions)(ii, 2) = (*ghost_positions)(i, 2) + offset[2];
          (*ghost_types)(ii) = typemap.at(types(i));
          ++ii;
        }
      }
}



unique_ptr<Box> Box::delete_atom(unsigned i, const string& name) const {
  if (i >= natoms_)
    throw runtime_error("Deleting non-existant atom.");
  unsigned new_natoms = natoms_ - 1;
  unique_ptr< Array2D<double> > new_positions =
    make_unique< Array2D<double> >(new_natoms, 3);
  unique_ptr< Array1DInit<string> > new_types =
    make_unique< Array1DInit<string> >(new_natoms);
  // Copy all except atom i.
  for (unsigned j = 0; j < i; ++j) {
    for (unsigned dim = 0; dim < 3; ++dim)
      (*new_positions)(j, dim) = positions(j, dim);
    (*new_types)(j) = types(j);
  }
  for (unsigned j = i+1; j < natoms_; ++j) {
    for (unsigned dim = 0; dim < 3; ++dim)
      (*new_positions)(j-1, dim) = positions(j, dim);
    (*new_types)(j-1) = types(j);
  }
  return make_unique<Box>(a, b, c, periodic[0], periodic[1], periodic[2],
                          move(new_positions), move(new_types),
                          name);
}



unique_ptr<Box> Box::repeat(unsigned repeat_a,
                            unsigned repeat_b,
                            unsigned repeat_c,
                            const std::string& name) const {
  unsigned new_natoms = natoms_ * repeat_a * repeat_b * repeat_c;
  unique_ptr< Array2D<double> > new_positions =
    make_unique< Array2D<double> >(new_natoms, 3);
  unique_ptr< Array1DInit<string> > new_types =
    make_unique< Array1DInit<string> >(new_natoms);
  unsigned new_idx = 0;
  for (unsigned i = 0; i < repeat_a; ++i) {
    for (unsigned j = 0; j < repeat_b; ++j) {
      for (unsigned k = 0; k < repeat_c; ++k) {
        Vec3D<double> offset(0.0, 0.0, 0.0);
        offset += double(i) * a;
        offset += double(j) * b;
        offset += double(k) * c;
        for (unsigned idx = 0; idx < natoms_; ++idx) {
          for (unsigned dim = 0; dim < 3; ++dim) {
            (*new_positions)(new_idx, dim) = positions(idx, dim) + offset[dim];
          }
          (*new_types)(new_idx) = types(idx);
          ++new_idx;
        }
      }
    }
  }
  return make_unique<Box>(a*double(repeat_a),
                          b*double(repeat_b),
                          c*double(repeat_c),
                          periodic[0], periodic[1], periodic[2],
                          move(new_positions), move(new_types),
                          name);
}


double Box::calc_dist(int i, int j, double& dx, double& dy, double& dz) const {
  if (!ghost_positions)
    throw runtime_error("Ghosts not initialized!");
  dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
  dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
  dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
  return sqrt(dx*dx + dy*dy + dz*dz);
}



void Box::write_to(ostream& output) const {
  // Number of atoms.
  output << natoms_ << "\n";
  // Extended XYZ uses the comment line for metadata like the box shape.
  output << "Lattice=\""
         << a[0] << " " << a[1] << " " << a[2] << " "
         << b[0] << " " << b[1] << " " << b[2] << " "
         << c[0] << " " << c[1] << " " << c[2] << "\" "
         << "Properties=species:S:1:pos:R:3\n";
  // The atoms.
  for (unsigned i = 0; i != natoms_; ++i)
    output << types(i) << " "
           << positions(i,0) << " "
           << positions(i,1) << " "
           << positions(i,2) << "\n";
  output << flush;
}


void Box::scale(double factor_a, double factor_b, double factor_c,
                const map<string,int>& typemap) {
  // Get transformation matrix.
  Array2D<double> T(3,3); T = 0.0;
  const double vol = calc_volume();
  const Vec3D<double> bxc = cross(b_, c_);
  const Vec3D<double> cxa = cross(c_, a_);
  const Vec3D<double> axb = cross(a_, b_);

  for (unsigned i = 0; i != 3; ++i)
    for (unsigned j = 0; j != 3; ++j)
      T(i,j) = ( a_[i]*factor_a*bxc[j] + b_[i]*factor_b*cxa[j] + c_[i]*factor_c*axb[j] ) / vol;

  // Update box.
  a_ *= factor_a;
  b_ *= factor_b;
  c_ *= factor_c;

  // Update positions
  for (unsigned atom = 0; atom != natoms_; ++atom) {
    Vec3D<double> pos(0.0, 0.0, 0.0);
    for (unsigned i = 0; i != 3; ++i)
      for (unsigned j = 0; j != 3; ++j) {
        pos[i] += T(i,j)*positions(atom,j);
      }
    // Assign back.
    for (unsigned dim = 0; dim != 3; ++dim)
      positions(atom,dim) = pos[dim];
  }

  // Update ghosts etc.
  update_ghosts(typemap);
}


void Box::deform(Voigt6<double> defmatrix,
                 const map<string,int>& typemap) {
  // Throw away very small shear.
  if (abs(defmatrix(3)) < DELTA)
    defmatrix(3) = 0.0;
  if (abs(defmatrix(4)) < DELTA)
    defmatrix(4) = 0.0;
  if (abs(defmatrix(5)) < DELTA)
    defmatrix(5) = 0.0;

  // Update box.
  Vec3D<double> new_a(0.0, 0.0, 0.0);
  Vec3D<double> new_b(0.0, 0.0, 0.0);
  Vec3D<double> new_c(0.0, 0.0, 0.0);
  for (unsigned i = 0; i != 3; ++i)
    for (unsigned j = 0; j != 3; ++j) {
      new_a[i] += defmatrix(i,j)*a_[j];
      new_b[i] += defmatrix(i,j)*b_[j];
      new_c[i] += defmatrix(i,j)*c_[j];
    }
  a_ = new_a; b_ = new_b; c_ = new_c;
  box_side_lengths_[0] = a_.abs();
  box_side_lengths_[1] = b_.abs();
  box_side_lengths_[2] = c_.abs();

  // Update positions.
  for (unsigned atom = 0; atom != natoms_; ++atom) {
    Vec3D<double> pos(0.0, 0.0, 0.0);
    for (unsigned i = 0; i != 3; ++i)
      for (unsigned j = 0; j != 3; ++j) {
        pos[i] += defmatrix(i,j)*positions(atom,j);
      }
    // Assign back.
    for (unsigned dim = 0; dim != 3; ++dim)
      positions(atom,dim) = pos[dim];
  }

  // Update ghosts etc.
  update_ghosts(typemap);
}



unsigned Box::atoms_per_unit_cell(const string& lattice, bool cubic) {

  if (lattice == "dimer")
    return 2;

  if (lattice == "sc")
    return 1;

  if (lattice == "bcc")
    return cubic ? 2 : 1;

  if (lattice == "fcc")
    return cubic ? 4 : 1;

  if (lattice == "diamond")
    return cubic ? 8 : 2;

  if (lattice == "B1")
    return cubic ? 8 : 2;

  if (lattice == "B2")
    return 2;

  if (lattice == "B3")
    return cubic ? 8 : 2;

  throw runtime_error("unknown lattice");
}

unsigned Box::species_per_unit_cell(const string& lattice) {
  if (lattice == "sc" || lattice == "bcc" || lattice == "fcc"
      || lattice == "diamond")
    return 1;
  else if (lattice == "dimer" || lattice == "B1" || lattice == "B2"
           || lattice == "B3")
    return 2;
  else
    throw runtime_error("unknown lattice");
}

Vec3D<unsigned> Box::calc_number_of_ghost_shells(double cutoff) const {
  const Vec3D<double> bc = cross(b_, c_);
  const Vec3D<double> ac = cross(a_, c_);
  const Vec3D<double> ab = cross(a_, b_);
  const double V = abs(dot(a_, bc));
  const double da = V / bc.abs();
  const double db = V / ac.abs();
  const double dc = V / ab.abs();
  return Vec3D<unsigned>(periodic_[0] ? ceil(cutoff/da) : 0,
                         periodic_[1] ? ceil(cutoff/db) : 0,
                         periodic_[2] ? ceil(cutoff/dc) : 0);
}





// Compute /////////////////////////////////////////////////////////////


Compute::Compute(unique_ptr<Box> box, const string& modelname,
                 KIM::LengthUnit length_unit,
                 KIM::EnergyUnit energy_unit,
                 KIM::ChargeUnit charge_unit,
                 KIM::TemperatureUnit temperature_unit,
                 KIM::TimeUnit time_unit)
  : box_(move(box)), modelname_(modelname)
{
  int error;

  // Create KIM model. /////////////////////////////////////////////////
  int requested_units_accepted;
  error = KIM::Model::Create(KIM::NUMBERING::zeroBased, // TODO: test this??    
                             length_unit,
                             energy_unit,
                             charge_unit,
                             temperature_unit,
                             time_unit,
                             modelname,
                             &requested_units_accepted,
                             &model);
  if (error)
    throw runtime_error("KIM model creation failed.");

  if (!requested_units_accepted)
    throw runtime_error("KIM model does not accept my units.");

  // Check routines and register support. //////////////////////////////
  int n_model_routines;
  KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(&n_model_routines);
  for (int i = 0; i < n_model_routines; ++i) {
    KIM::ModelRoutineName model_routine_name;
    KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &model_routine_name);

    int present, required;
    error = model->IsRoutinePresent(model_routine_name, &present, &required);
    if (error)
      throw runtime_error("KIM model->IsRoutinePresent() failed for "
                          + model_routine_name.ToString() + ".");

    if (model_routine_name
        == KIM::MODEL_ROUTINE_NAME::Create
        ||
        model_routine_name
        == KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate
        ||
        model_routine_name
        == KIM::MODEL_ROUTINE_NAME::Compute
        ||
        model_routine_name
        == KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy
        ||
        model_routine_name
        == KIM::MODEL_ROUTINE_NAME::Destroy) {
      // Those are required by the API, so let's just check they are there.
      if (!present)
        throw runtime_error("The required routine \""
                            + model_routine_name.ToString()
                            + "\" is missing.");
    } else if (model_routine_name
               == KIM::MODEL_ROUTINE_NAME::Extension) {
      // Not supported by me.
      if (present && required)
        throw runtime_error("The model requires an extension. "
                            "This is not supported.");
    } else if (model_routine_name
               == KIM::MODEL_ROUTINE_NAME::Refresh) {
      has_reinit = present;
    } else if (model_routine_name
               == KIM::MODEL_ROUTINE_NAME::WriteParameterizedModel) {
      has_write_params = present;
    } else {
      // Unknown.
      if (present && required)
        throw runtime_error("The unknown routine \""
                            + model_routine_name.ToString()
                            + "\" is required by the model.");
    }
  }

  // TODO: We could get back the units and check, but do we need to?      

  // Check species. ////////////////////////////////////////////////////
  for (int i = 0; i < box_->types.extent(0); ++i) {
    int species_supported;
    int species_code;
    // Get code from model.
    const KIM::SpeciesName sn(box_->types(i));
    error = model->GetSpeciesSupportAndCode(sn, &species_supported,
                                            &species_code);
    if (error || !species_supported)
      throw runtime_error("Species \""
                          + box_->types(i)
                          + "\" from box is not supported by the model.");
    partcl_type_codes[box_->types(i)] = species_code;
    partcl_type_names[species_code] = box_->types(i);
  }

  // Init KIM compute arguments. ///////////////////////////////////////
  error = model->ComputeArgumentsCreate(&compute_arguments);

  if (error)
    throw runtime_error("KIM model->ComputeArgumentsCreate() failed.");

  // Check and set the compute arguments.
  int n_comp_args;
  KIM::COMPUTE_ARGUMENT_NAME::GetNumberOfComputeArgumentNames(&n_comp_args);
  for (int i = 0; i < n_comp_args; ++i) {
    KIM::ComputeArgumentName compute_argument_name;
    KIM::COMPUTE_ARGUMENT_NAME::GetComputeArgumentName(i,
                                                       &compute_argument_name);
    // TODO: get its type and support/required status and check if we are compatible     

    // Set the arguments.
    if (compute_argument_name
        == KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles) {
      error = compute_arguments->SetArgumentPointer(compute_argument_name,
                                                    &nall_kim);
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::particleContributing) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::coordinates) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::partialEnergy) {
      error = compute_arguments->SetArgumentPointer(compute_argument_name,
                                                    &energy);
      has_energy = true; // TODO: get from API    
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
      has_particleEnergy = true; // TODO: get from API    
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::partialForces) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
      has_forces = true; // TODO: get from API    
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::partialVirial) {
      error = compute_arguments->SetArgumentPointer(compute_argument_name,
                                                    &virial(0));
      has_virial = true; // TODO: get from API    
    } else if (compute_argument_name
               == KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial) {
      // The length of this array is only determined after neighbor
      // lists have been calculated and the argument will thus be set
      // by update_neighbor_list().
      has_particleVirial = true; // TODO: get from API    
    } else {
      // TODO: complain if something else required    
    }
    if (error)
      throw runtime_error("KIM model->SetArgumentPointer() failed for "
                          + compute_argument_name.ToString() + ".");
  }

  // Check and set the callbacks.
  int n_comp_callb;
  KIM::COMPUTE_CALLBACK_NAME::GetNumberOfComputeCallbackNames(&n_comp_callb);
  for (int i = 0; i < n_comp_callb; ++i) {
    KIM::ComputeCallbackName compute_callback_name;
    KIM::COMPUTE_CALLBACK_NAME::GetComputeCallbackName(i,
                                                       &compute_callback_name);
    // TODO: get its support/required status and check if we are compatible     

    if (compute_callback_name
        == KIM::COMPUTE_CALLBACK_NAME::GetNeighborList) {
      error = compute_arguments->SetCallbackPointer(compute_callback_name,
                                                    KIM::LANGUAGE_NAME::cpp,
                                                    (KIM::Function*) &get_neigh,
                                                    this);
    } else if (compute_callback_name
               == KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm) {
      error =
        compute_arguments->SetCallbackPointer(compute_callback_name,
                                              KIM::LANGUAGE_NAME::cpp,
                                              (KIM::Function*) &process_dEdr,
                                              this);
      has_process_dEdr = true;
    } else {
      //TODO complain if required  
    }
    if (error)
      throw runtime_error("KIM model->SetCallbackPointer() failed for "
                          + compute_callback_name.ToString() + ".");
  }

  // Get cutoff(s) and create neighbor lists. //////////////////////////
  model->GetInfluenceDistance(&cutoff);

  //TODO: implement that shit, for now we just have one neighlist :-)     
  /*
  model->GetNeighborListPointers(&NUMBER_OF_NEIGHBORLISTS,   
                                 &CUTOFFS,     
                                 &WILL_REQUEST_NEIGHBORS_OF_GHOSTS);    
  */

  // Now that we know the cutoff, we calculate neighbor lists and
  // ghost atoms.  Then we allocate memory for variable length data
  // and pass it to KIM.
  update_neighbor_list();

  // Get parameters of the model and store 'em. ////////////////////////

  if (has_reinit) {
    int n_params;
    model->GetNumberOfParameters(&n_params);

    for (int i; i < n_params; ++i) {
      KIM::DataType data_type;
      const string * param_name;
      const string * param_desc;
      int size;
      model->GetParameterMetadata(i, &data_type, &size,
                                  &param_name, &param_desc);
      // Fancy way to do free_parameter_map[param_name] = fp; but
      // without assignment or copy.
      free_parameter_map.emplace(piecewise_construct,
                                 forward_as_tuple(*param_name),
                                 forward_as_tuple(*param_name,
                                                  *param_desc,
                                                  data_type,
                                                  size,
                                                  i));
    }
  }
}

Compute::~Compute() {
  KIM::Model::Destroy(&model);
}


void Compute::compute() {
  // Reset arrays to zero. Dunno if needed, but with ghost atoms I am
  // unsure.
  /*
  if (has_forces)
    for (unsigned i = 0; i < forces.size(); ++i)
      forces[i] = 0;
  if (has_particleEnergy)
    for (unsigned i = 0; i < particleEnergy.size(); ++i)
      particleEnergy[i] = 0;
  if (has_particleVirial)
    for (unsigned i = 0; i < particleVirial.size(); ++i)
      particleVirial[i] = 0;
  */
  for (unsigned dim = 0; dim < 6; ++dim)
    virial_from_dEdr(dim) = 0.0;
  // Compute.
  const int error = model->Compute(compute_arguments);
  if (error)
    throw runtime_error(string("Error in KIM's compute function in line ")
                        + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Before reducing the ghost atoms, calculate the global virial from
  // forces. This can be used to verify the global virial returned by
  // the model.
  if (has_forces) {
    global_virial_from_forces_xx = 0.0;
    global_virial_from_forces_yy = 0.0;
    global_virial_from_forces_zz = 0.0;
    global_virial_from_forces_yz = 0.0;
    global_virial_from_forces_xz = 0.0;
    global_virial_from_forces_xy = 0.0;
    const double* ghost_pos = box_->get_positions_ptr();
    for (unsigned i = 0; i < box_->nall; ++i) {
      global_virial_from_forces_xx -= ghost_pos[3*i + 0] * forces[3*i + 0];
      global_virial_from_forces_yy -= ghost_pos[3*i + 1] * forces[3*i + 1];
      global_virial_from_forces_zz -= ghost_pos[3*i + 2] * forces[3*i + 2];
      global_virial_from_forces_yz -= ghost_pos[3*i + 1] * forces[3*i + 2];
      global_virial_from_forces_xz -= ghost_pos[3*i + 0] * forces[3*i + 2];
      global_virial_from_forces_xy -= ghost_pos[3*i + 0] * forces[3*i + 1];
    }
  }

  // Tricks to fix values when ghost atoms are in the game.
  if (box_->nghosts)
    for (unsigned i = box_->natoms; i < box_->nall; ++i) {
      const unsigned central = i % box_->natoms;
      if (has_particleEnergy) {
        particleEnergy[central] += particleEnergy[i];
      }
      if (has_forces) {
        forces[3*central + 0] += forces[3*i + 0];
        forces[3*central + 1] += forces[3*i + 1];
        forces[3*central + 2] += forces[3*i + 2];
      }
      if (has_particleVirial) {
        particleVirial[6*central + 0] += particleVirial[6*i + 0];
        particleVirial[6*central + 1] += particleVirial[6*i + 1];
        particleVirial[6*central + 2] += particleVirial[6*i + 2];
        particleVirial[6*central + 3] += particleVirial[6*i + 3];
        particleVirial[6*central + 4] += particleVirial[6*i + 4];
        particleVirial[6*central + 5] += particleVirial[6*i + 5];
      }
    }
}


int Compute::get_neigh(void * const compute_ptr,
                       const int number_of_neighlists,
                       const double * const cutoffs,
                       const int neighbor_list_idx,
                       const int i,
                       int * const n_neighs,
                       const int ** const neighlist) {
  // TODO: implement different neighbor lists if requested        

  Compute& c = *( (Compute *)compute_ptr );

  if (i < 0 || i >= c.nall_kim) {
    // Out of range access.
    return 1;
  } else if (i < c.nall_kim) {
    // Not a ghost atom.
    const vector<int>& neighbors = c.box_->get_neighbors(i);
    *n_neighs = neighbors.size();
    *neighlist = &neighbors[0];
  } else {
    // Ghost atom. Has no neighbors.
    *n_neighs = 0;
    *neighlist = nullptr;
  }

  return 0;
}


int Compute::process_dEdr(void * const compute_ptr,
                          const double dEdr,
                          const double r,
                          const double * const vec_r,
                          const int i,
                          const int j) {
  Compute& c = *( (Compute *)compute_ptr );

  // Compute global virial from dEdr.
  const double v = dEdr / r;
  for (unsigned ii = 0; ii < 3; ++ii)
    for (unsigned jj = 0; jj < 3; ++jj)
      c.virial_from_dEdr(ii,jj) +=
        ((ii == jj) ? 1.0 : 0.5) * v * vec_r[ii] * vec_r[jj];

  return 0;
}


void Compute::update_neighbor_list(bool force_ptr_update) {
  // TODO: For now the skin is 0.
  const bool arrays_changed = box_->update_neighbor_list(cutoff, 0.0,
                                                         partcl_type_codes);
  if (force_ptr_update || arrays_changed) {
    const unsigned n = box_->nall;
    nall_kim = n; // update the number of atoms variable exposed to KIM.
    if (has_forces) forces.resize(3 * n);
    if (has_particleEnergy) particleEnergy.resize(n);
    if (has_particleVirial) particleVirial.resize(6 * n);
    int error =
      compute_arguments->SetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
        box_->get_types_ptr())
      ||
      compute_arguments->SetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
        box_->get_contributing_ptr())
      ||
      compute_arguments->SetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::coordinates,
        box_->get_positions_ptr());
    if (has_particleEnergy)
      error = error || compute_arguments->SetArgumentPointer(
                         KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
                         &particleEnergy[0]);
    if (has_forces)
      error = error || compute_arguments->SetArgumentPointer(
                         KIM::COMPUTE_ARGUMENT_NAME::partialForces,
                         &forces[0]);
    if (has_particleVirial)
      error = error || compute_arguments->SetArgumentPointer(
                         KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
                         &particleVirial[0]);
    if (error)
      throw runtime_error(string("Error in KIM's SetArgumentPointer in line ")
                          + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
  }
}

template<typename T>
void Compute::set_parameter_impl(const FreeParam& param,
                                 const unsigned index,
                                 T value, bool reinit) {
  // Check if the model supports reinit.
  if (!has_reinit)
    throw runtime_error("model does not support changing parameters.");
  // Check parameter type.
  if (param.data_type == KIM::DATA_TYPE::Integer
      &&
      !is_same<T, int>::value)
    throw runtime_error("parameter \"" + param.name + "\" is an int!");
  else if (param.data_type == KIM::DATA_TYPE::Double
           &&
           !is_same<T, double>::value)
    throw runtime_error("parameter \"" + param.name + "\" is a double!");
  // Check index.
  if (index >= param.size)
    throw runtime_error("index " + to_string(index) + " out of bounds [0; "
                        + to_string(param.size-1) + "].");
  // Set it.
  int error = model->SetParameter(param.kim_index, index, value);
  if (error)
    throw runtime_error("KIM model SetParameter() returned an error");
  // Update KIM model if requested.
  if (reinit) {
    error = model->ClearThenRefresh();
    if (error)
      throw runtime_error("KIM model ClearThenRefresh() returned an error");
    // Get cutoff.
    model->GetInfluenceDistance(&cutoff);
    //TODO: if multiple neighbor lists supported, updated them, too    
  }
}

void Compute::set_parameter(const string& param_name,
                            const unsigned index,
                            double value, bool reinit) {
  set_parameter_impl<double>(get_param_obj(param_name), index, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const unsigned index,
                            int value, bool reinit) {
  set_parameter_impl<int>(get_param_obj(param_name), index, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const string& species1,
                            const string& species2,
                            double value, bool reinit) {
  const FreeParam& param = get_param_obj(param_name);
  const unsigned index = conv_index(get_particle_type_code(species1),
                                    get_particle_type_code(species2),
                                    param.size);
  set_parameter_impl<double>(param, index, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const string& species1,
                            const string& species2,
                            int value, bool reinit) {
  const FreeParam& param = get_param_obj(param_name);
  const unsigned index = conv_index(get_particle_type_code(species1),
                                    get_particle_type_code(species2),
                                    param.size);
  set_parameter_impl<int>(param, index, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const string& species1,
                            const string& species2,
                            const string& species3,
                            double value, bool reinit) {
  const FreeParam& param = get_param_obj(param_name);
  const unsigned index = conv_index(get_particle_type_code(species1),
                                    get_particle_type_code(species2),
                                    get_particle_type_code(species3),
                                    param.size);
  set_parameter_impl<double>(param, index, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const string& species1,
                            const string& species2,
                            const string& species3,
                            int value, bool reinit) {
  const FreeParam& param = get_param_obj(param_name);
  const unsigned index = conv_index(get_particle_type_code(species1),
                                    get_particle_type_code(species2),
                                    get_particle_type_code(species3),
                                    param.size);
  set_parameter_impl<int>(param, index, value, reinit);
}



template<typename T>
T Compute::get_parameter_impl(const string& param_name,
                              const unsigned index) {
  // Get parameter data.
  const auto it = free_parameter_map.find(param_name);
  if (it == free_parameter_map.end())
    throw runtime_error("Unknown free parameter: " + param_name);
  const FreeParam& param = it->second;
  // Check parameter type.
  if (param.data_type == KIM::DATA_TYPE::Integer
      &&
      !is_same<T, int>::value)
    throw runtime_error("parameter \"" + param_name + "\" is an int!");
  else if (param.data_type == KIM::DATA_TYPE::Double
           &&
           !is_same<T, double>::value)
    throw runtime_error("parameter \"" + param_name + "\" is a double!");
  // Check index.
  if (index >= param.size)
    throw runtime_error("index " + to_string(index) + " out of bounds [0; "
                        + to_string(param.size-1) + "].");
  // Get it.
  T retval;
  const int error = model->GetParameter(param.kim_index, index, &retval);
  if (error)
    throw runtime_error("KIM model GetParameter() returned an error");
  return retval;
}

int Compute::get_parameter_int(const string& param_name,
                               const unsigned index) {
  return get_parameter_impl<int>(param_name, index);
}

double Compute::get_parameter_double(const string& param_name,
                                     const unsigned index) {
  return get_parameter_impl<double>(param_name, index);
}


// Fitting /////////////////////////

double Compute::obj_func_box(const vector<double>& x, vector<double>& grad,
                             void* f_data) {
  Compute& c = *(static_cast<Compute*>(f_data));
  if (x.size() != 3 && x.size() != 1)
    throw runtime_error("This objective function only works with "
                        "1 or 3 parameters.");
  // Scale box, we assume the neighbor list doesn't change (TODO?).
  if (x.size() == 3)
    c.box_->scale_to(x[0], x[1], x[2], c.partcl_type_codes);
  else {
    const double factor = x[0] / c.box_->box_side_lengths[0];
    c.box_->scale(factor, c.partcl_type_codes);
  }
  c.compute();
  ++c.fit_counter;
  // Gradient.
  if (!grad.empty())
    throw runtime_error("Gradient not supported for box optimization.");
  // Objective function value is energy.
  return c.energy;
}

double Compute::optimize_box(double ftol_abs, unsigned maxeval,
                             bool isotropic) {
  // Check if optimizing box is supported. Box optimization needs
  // energy, so check if that is provided.
  if (!has_energy)
    throw runtime_error("The model does not provide energy calculation "
                        "which is needed to optimize the box.");
  //
  vector<double> lengths;
  lengths.push_back(box_->box_side_lengths[0]);
  if (!isotropic) {
    lengths.push_back(box_->box_side_lengths[1]);
    lengths.push_back(box_->box_side_lengths[2]);
  }
  static const vector<double> lb1 = { 0.0 }; // Lengths may not be negative.
  static const vector<double> lb3 = { 0.0, 0.0, 0.0 };
  // TODO: We use a gradient-free algorithm for now, although the
  // gradient should be obtainable from the virial (is it even exactly
  // the virial?).
  double obj_val;
  nlopt::opt optimizer(nlopt::LN_SBPLX, isotropic ? 1 : 3);
  optimizer.set_min_objective(Compute::obj_func_box, this);
  optimizer.set_lower_bounds(isotropic ? lb1 : lb3);
  optimizer.set_initial_step(0.05); // This may not be too big!  TODO:
                                    // user-definable or better
                                    // heuristics?
  fit_counter = 0;
  optimizer.set_maxeval(maxeval);
  optimizer.set_ftol_abs(ftol_abs);
  optimizer.optimize(lengths, obj_val);
  return obj_val;
}


double Compute::obj_func_pos(const vector<double>& x, vector<double>& grad,
                             void* f_data) {
  Compute& c = *(static_cast<Compute*>(f_data));
  // Write back all positions.
  double* pos = &(c.box_->positions(0,0));
  for (unsigned i = 0; i != c.box_->natoms; ++i)
    pos[i] = x[i];
  c.box_->update_ghosts(c.partcl_type_codes);
  c.compute();
  ++c.fit_counter;
  // Fill gradient.
  if (!grad.empty()) {
    const double* f = &c.forces[0];
    for (unsigned i = 0; i != c.box_->natoms; ++i)
      grad[i] = -f[i];
  }
  // Objective function value is energy.
  return c.energy;
}

double Compute::optimize_positions(double ftol_abs, unsigned maxeval) {
  // Check if optimizing box is supported. Box optimization needs
  // energy, so check if that is provided.
  if (!has_energy)
    throw runtime_error("The model does not provide energy calculation "
                        "which is needed to optimize atom positions.");
  // Get coordinates as vector.
  Array2D<double>& pos = box_->positions;
  const int nparams = pos.extent(0)*pos.extent(1);
  const double* p = &pos(0,0);
  vector<double> pos_vec(p, p+nparams);
  // Other good algorithms:
  //   * NLOPT_LD_LBFGS
  //   * NLOPT_LD_TNEWTON_PRECOND_RESTART
  //   * NLOPT_LD_VAR2
  nlopt::algorithm algo;
  if (has_forces)
    // We can use an algorithm with gradients.
    algo = nlopt::LD_TNEWTON_PRECOND_RESTART;
  else {
    // We must use a gradient-free algorithm.
    cerr << "WARNING: model does not compute forces, "
         << "atomic position optimization will be slow." << endl;
    algo = nlopt::LN_SBPLX;
  }
  double obj_val;
  nlopt::opt optimizer(algo, nparams);
  optimizer.set_min_objective(Compute::obj_func_pos, this);
  fit_counter = 0;
  optimizer.set_maxeval(maxeval);
  optimizer.set_ftol_abs(ftol_abs);
  optimizer.optimize(pos_vec, obj_val);
  return obj_val;
}



void Compute::switch_boxes(Compute& other) {
  // Check model.
  /* I don't think I need this; the supported species are checked and that's it
  if (modelname_ != other.modelname_)
    throw runtime_error("KIM models do not match.");
  */
  // Swap boxes.
  box_.swap(other.box_);
  // Reset KIM data on this object.
  update_kim_after_box_change();
  // Reset KIM data on other object.
  other.update_kim_after_box_change();
}

unique_ptr<Box> Compute::change_box(unique_ptr<Box> new_box) {
  box_.swap(new_box);
  update_kim_after_box_change();
  return new_box; // Which is now the old box.
}


void Compute::update_kim_after_box_change() {
  // Check and update species.
  partcl_type_codes.clear();
  partcl_type_names.clear();
  for (int i = 0; i < box_->types.extent(0); ++i) {
    int species_supported;
    int species_code;
    // Get code from model.
    const KIM::SpeciesName sn(box_->types(i));
    const int error =
      model->GetSpeciesSupportAndCode(sn, &species_supported, &species_code);
    if (error || !species_supported)
      throw runtime_error("Species \""
                          + box_->types(i)
                          + "\" from box is not supported by the model.");
    partcl_type_codes[box_->types(i)] = species_code;
    partcl_type_names[species_code] = box_->types(i);
  }
  // Also force updating all pointers and resizing all output array since
  // the box changed.
  update_neighbor_list(true);
}
