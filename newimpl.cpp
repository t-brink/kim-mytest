#include <iostream>
#include <cstdio>
#include <utility>
#include <cmath>

#include "mytest.hpp"

using namespace std;
using namespace mytest;

Box::Box(const Vec3D<double>& a_in,
         const Vec3D<double>& b_in,
         const Vec3D<double>& c_in,
         bool periodic_a, bool periodic_b, bool periodic_c,
         unique_ptr< Array2D<double> > coordinates,
         unique_ptr< Array1D<int> > types_in,
         bool neigh_list, const string& name)
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
    neigh_list_(neigh_list ? natoms_ : 0),
    neigh_rvec_(neigh_list ? natoms_ : 0)
{
  if (positions.extent(0) != types.extent(0))
    throw runtime_error("types must have the same length as the "
                        "first dimension of positions.");
  if (positions.extent(1) != 3)
    throw runtime_error("positions must be a n*3 array.");

  // TODO: neigh_list, ghosts?       
}

Box::Box(const std::string& lattice, double lattice_const, bool cubic,
         unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
         bool periodic_a, bool periodic_b, bool periodic_c,
         const std::vector<int>& types_in,
         bool neigh_list, const std::string& name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(atoms_per_unit_cell(lattice, cubic) // Public data.
              * repeat_a * repeat_b * repeat_c, 3),
    types(positions.extent(0)),
    periodic_(periodic_a, periodic_b, periodic_c), // Private data.
    natoms_(positions.extent(0)), nghosts_(0), nall_(natoms_),
    name_(name)
{
  if (types_in.size() != species_per_unit_cell(lattice))
    throw runtime_error("Wrong number of types for this lattice.");

  // The unit cell vectors.
  Vec3D<double> unit_a;
  Vec3D<double> unit_b;
  Vec3D<double> unit_c;
  typedef pair<Vec3D<double>,int> t_atom;
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
    throw runtime_error("not implemented.");
  } else if (lattice == "B1") {
    throw runtime_error("not implemented.");
  } else if (lattice == "B2") {
    throw runtime_error("not implemented.");
  } else if (lattice == "B3") {
    throw runtime_error("not implemented.");
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
  if (neigh_list) {
    neigh_list_.resize(natoms_);
    neigh_rvec_.resize(natoms_);
  }

  // TODO: neigh_list, ghosts?       
}

bool Box::update_neighbor_list(double cutoff, double skin) {
  const double cut = (1 + skin) * cutoff;
  ghost_shells = calc_number_of_ghost_shells(cut);
  const unsigned new_nghosts = natoms_ * ((2*ghost_shells[0]+1)
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
    reallocated = true;
  }

  // Update ghosts.
  update_ghosts();

  // Clear neighbor lists.
  for (auto& l : neigh_list_)
    l.clear();
  for (auto& l : neigh_rvec_)
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
        neigh_rvec_[i].push_back(dx);
        neigh_rvec_[i].push_back(dy);
        neigh_rvec_[i].push_back(dz);
        neigh_rvec_[j].push_back(-dx);
        neigh_rvec_[j].push_back(-dy);
        neigh_rvec_[j].push_back(-dz);
      }
    }
    // Iterate over ghosts.
    for (unsigned j = natoms_; j != nall_; ++j) {
      const double dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
      const double dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
      const double dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
      if (dx*dx + dy*dy + dz*dz < cutsq) {
        neigh_list_[i].push_back(j);
        neigh_rvec_[i].push_back(dx);
        neigh_rvec_[i].push_back(dy);
        neigh_rvec_[i].push_back(dz);
      }
    }
  }

  return reallocated;
}

void Box::update_ghosts() {
  for (unsigned i = 0; i != natoms_; ++i) {
    (*ghost_positions)(i, 0) = positions(i, 0);
    (*ghost_positions)(i, 1) = positions(i, 1);
    (*ghost_positions)(i, 2) = positions(i, 2);
    (*ghost_types)(i) = types(i);
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
          (*ghost_positions)(ii, 0) = positions(i, 0) + offset[0];
          (*ghost_positions)(ii, 1) = positions(i, 1) + offset[1];
          (*ghost_positions)(ii, 2) = positions(i, 2) + offset[2];
          (*ghost_types)(ii) = types(i);
          ++ii;
        }
      }
}


double Box::calc_dist(int i, int j, double& dx, double& dy, double& dz) const {
  if (!ghost_positions)
    throw runtime_error("Ghosts no initialized!");
  dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
  dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
  dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
  return sqrt(dx*dx + dy*dy + dz*dz);
}



void Box::write_to(ostream& output) const {
  // Comment line
  output << name_ << "\n";
  // Scaling factor, always 1.
  output << 1.0 << "\n";
  // Box vectors;
  output << a[0] << " " << a[1] << " " << a[2] << "\n";
  output << b[0] << " " << b[1] << " " << b[2] << "\n";
  output << c[0] << " " << c[1] << " " << c[2] << "\n";
  // Types. TODO: for now only one type.
  output << natoms_ << "\n";
  // Cartesian coordinates.
  output << "cartesian\n";
  for (unsigned i = 0; i != natoms_; ++i)
    output << positions(i,0) << " "
           << positions(i,1) << " "
           << positions(i,2) << "\n";
  output << endl;
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

unsigned Box::species_per_unit_cell(const std::string& lattice) {
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


int main() {
  vector<int> types;
  types.push_back(1);
  //types.push_back(1);
  Box b("fcc", 3.940, true, 10, 10, 10, false, false, false,
        types, true, "box");

  b.update_neighbor_list(2.92, 0.0);

  const int i = 0;
  cout << "Neighbors of atom " << i << ":";
  for (const auto& j : b.get_neighbors(i)) {
    printf("  %6d", j);
  }
  cout << endl;
  cout << "    distance from " << i << ":";
  for (const auto& j : b.get_neighbors(i)) {
    printf("  %6.3f", b.calc_dist(i,j));
  }
  cout << endl;

  b.write_to("dump");

  return 0;
}
