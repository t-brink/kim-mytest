//#include <iostream>
//#include <cstdio>
#include <utility>
#include <cmath>
#include <chrono>

#include <KIM_API_C.h>
#include <KIM_API_status.h>
#include <nlopt.hpp>

#include "mytestcore.hpp"
#include "utils.hpp"

using namespace std;
using namespace mytest;

const double DELTA = 1e-10; // For double comparison.

// Box /////////////////////////////////////////////////////////////////

Box::Box(const Vec3D<double>& a_in,
         const Vec3D<double>& b_in,
         const Vec3D<double>& c_in,
         bool periodic_a, bool periodic_b, bool periodic_c,
         unique_ptr< Array2D<double> > coordinates,
         unique_ptr< Array1D<int> > types_in,
         KIMNeigh neighmode, const string& name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    kim_neighbor_mode(neighmode),
    natoms(natoms_), nghosts(nghosts_), nall(nall_),
    positions(move(*coordinates)), // Public data.
    types(move(*types_in)),
    box_side_lengths_(a_in.abs(), b_in.abs(), c_in.abs()), // Private data.
    a_(a_in), b_(b_in), c_(c_in),
    periodic_(periodic_a, periodic_b, periodic_c),
    natoms_(coordinates->extent(0)), nghosts_(0), nall_(natoms_),
    name_(name),
    neigh_list_(neighmode != KIM_cluster ? natoms_ : 0),
    neigh_rvec_(neighmode == KIM_neigh_rvec_f ? natoms_ : 0),
    neigh_rvec_shell_(neighmode == KIM_neigh_rvec_f ? natoms_ : 0)
{
  if (positions.extent(0) != types.extent(0))
    throw runtime_error("types must have the same length as the "
                        "first dimension of positions.");
  if (positions.extent(1) != 3)
    throw runtime_error("positions must be a n*3 array.");

  if (neighmode == KIM_mi_opbc_f
      && (!periodic_[0] || !periodic_[1] || !periodic_[2]
          || abs(a_in[1]) > DELTA || abs(a_in[2]) > DELTA
          || abs(b_in[0]) > DELTA || abs(b_in[2]) > DELTA
          || abs(c_in[0]) > DELTA || abs(c_in[1]) > DELTA))
      throw runtime_error("for MI_OPBC_F all dimensions must be periodic "
                          "and the box must be orthorhombic.");
}

Box::Box(const std::string& lattice, double lattice_const, bool cubic,
         unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
         bool periodic_a, bool periodic_b, bool periodic_c,
         const std::vector<int>& types_in,
         KIMNeigh neighmode, const std::string& name)
  : box_side_lengths(box_side_lengths_), // Public const references.
    a(a_), b(b_), c(c_),
    periodic(periodic_),
    kim_neighbor_mode(neighmode),
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
  if (neighmode == KIM_mi_opbc_f
      && (!periodic_[0] || !periodic_[1] || !periodic_[2]
          || abs(a_[1]) > DELTA || abs(a_[2]) > DELTA
          || abs(b_[0]) > DELTA || abs(b_[2]) > DELTA
          || abs(c_[0]) > DELTA || abs(c_[1]) > DELTA))
      throw runtime_error("for MI_OPBC_F all dimensions must be periodic "
                          "and the box must be orthorhombic.");
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
  if (neighmode != KIM_cluster)
    neigh_list_.resize(natoms_);
  if (neighmode == KIM_neigh_rvec_f) {
    neigh_rvec_.resize(natoms_);
    neigh_rvec_shell_.resize(natoms_);
  }
}

bool Box::update_neighbor_list(double cutoff, double skin) {
  const double cut = (1 + skin) * cutoff;
  if (kim_neighbor_mode == KIM_mi_opbc_f
      && (cut > 0.5 * box_side_lengths_[0]
          || cut > 0.5 * box_side_lengths_[1]
          || cut > 0.5 * box_side_lengths_[2]))
    throw runtime_error("Cutoff too big. When using MI_OPBC_F the cutoff "
                        "must be smaller than half the smallest box side "
                        "length.");
  ghost_shells = calc_number_of_ghost_shells(cut);
  Vec3D<unsigned> tmp_shells = ghost_shells;
  unsigned new_nghosts;
  switch (kim_neighbor_mode) {
  case KIM_cluster:
  case KIM_mi_opbc_f:
  case KIM_neigh_pure_f:
    // For these we store any ghost atoms that we may need.
    new_nghosts = natoms_ * ((2*ghost_shells[0]+1)
                             * (2*ghost_shells[1]+1)
                             * (2*ghost_shells[2]+1) - 1);
    break;
  case KIM_neigh_rvec_f:
    // Here the ghost atoms are only needed during neighbor list
    // calculation and are not stored.
    new_nghosts = 0;
    ghost_shells[0] = 0; ghost_shells[1] = 0; ghost_shells[2] = 0;
    break;
  default:
    throw runtime_error("unsupported neighbor list mode.");
  }
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
  switch(kim_neighbor_mode) {
  case KIM_cluster:
    break;
  case KIM_mi_opbc_f:
    {
      const double cutsq = cut * cut;
      for (unsigned i = 0; i != natoms_; ++i) {
        // Go over neighbors in the central cell.
        for (unsigned j = i+1; j != natoms_; ++j) {
          // Get distance vector.
          double dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
          if (abs(dx) > 0.5 * box_side_lengths_[0])
            dx -= (dx / abs(dx)) * box_side_lengths_[0];
          double dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
          if (abs(dy) > 0.5 * box_side_lengths_[1])
            dy -= (dy / abs(dy)) * box_side_lengths_[1];
          double dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
          if (abs(dz) > 0.5 * box_side_lengths_[2])
            dz -= (dz / abs(dz)) * box_side_lengths_[2];
          // Check cutoff.
          if (dx*dx + dy*dy + dz*dz < cutsq) {
            neigh_list_[i].push_back(j);
            neigh_list_[j].push_back(i);
          }
        }
      }
    }
    break;
  case KIM_neigh_pure_f:
    {
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
    }
    break;
  case KIM_neigh_rvec_f:
    {
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
            neigh_rvec_shell_[i].push_back(Vec3D<int>(0,0,0));
            neigh_rvec_shell_[j].push_back(Vec3D<int>(0,0,0));
          }
        }
        // Iterate over ghosts.
        const int alo = -static_cast<int>(tmp_shells[0]);
        const int ahi = tmp_shells[0];
        const int blo = -static_cast<int>(tmp_shells[1]);
        const int bhi = tmp_shells[1];
        const int clo = -static_cast<int>(tmp_shells[2]);
        const int chi = tmp_shells[2];
        for (int aa = alo; aa <= ahi; ++aa)
          for (int bb = blo; bb <= bhi; ++bb)
            for (int cc = clo; cc <= chi; ++cc) {
              if (aa == 0 && bb == 0 && cc == 0) continue;
              const Vec3D<double> offset =
                double(aa)*a + double(bb)*b + double(cc)*c;
              for (unsigned j = 0; j != natoms_; ++j) {
                const double dx =
                  (*ghost_positions)(j,0) - (*ghost_positions)(i,0) + offset[0];
                const double dy =
                  (*ghost_positions)(j,1) - (*ghost_positions)(i,1) + offset[1];
                const double dz =
                  (*ghost_positions)(j,2) - (*ghost_positions)(i,2) + offset[2];
                if (dx*dx + dy*dy + dz*dz < cutsq) {
                  neigh_list_[i].push_back(j);
                  neigh_rvec_[i].push_back(dx);
                  neigh_rvec_[i].push_back(dy);
                  neigh_rvec_[i].push_back(dz);
                  neigh_rvec_shell_[i].push_back(Vec3D<int>(aa,bb,cc));
                }
              }
            }
      }
    }
    break;
  default:
    throw runtime_error("unsupported neighbor list mode.");
  }

  return reallocated;
}

void Box::update_ghosts() {
  // Copy original atoms, enforcing periodic boundaries if the
  // neighbor mode is MI_OPBC_F.
  if (kim_neighbor_mode == KIM_mi_opbc_f)
    // Box is orthorhombic and periodic in all directions.
    for (unsigned i = 0; i != natoms_; ++i) {
      (*ghost_positions)(i, 0) = pmod(positions(i, 0), box_side_lengths_[0]);
      (*ghost_positions)(i, 1) = pmod(positions(i, 1), box_side_lengths_[1]);
      (*ghost_positions)(i, 2) = pmod(positions(i, 2), box_side_lengths_[2]);
      (*ghost_types)(i) = types(i);
    }
  else
    // Just assign without wrapping.
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
          (*ghost_positions)(ii, 0) = (*ghost_positions)(i, 0) + offset[0];
          (*ghost_positions)(ii, 1) = (*ghost_positions)(i, 1) + offset[1];
          (*ghost_positions)(ii, 2) = (*ghost_positions)(i, 2) + offset[2];
          (*ghost_types)(ii) = types(i);
          ++ii;
        }
      }
}


void Box::update_ghost_rvecs() {
  update_ghosts();
  if (kim_neighbor_mode == KIM_neigh_rvec_f) {
    for (unsigned i = 0; i != natoms_; ++i) {
      const unsigned jjmax = neigh_list_[i].size();
      const auto& neigh_i = neigh_list_[i];
      auto& rvec_i = neigh_rvec_[i];
      const auto& shells_i = neigh_rvec_shell_[i];
      for (unsigned jj = 0; jj != jjmax; ++jj) {
        const unsigned j = neigh_i[jj];
        const Vec3D<double> offset(shells_i[jj][0] * a[0]
                                   + shells_i[jj][1] * b[0]
                                   + shells_i[jj][2] * c[0],
                                   shells_i[jj][0] * a[1]
                                   + shells_i[jj][1] * b[1]
                                   + shells_i[jj][2] * c[1],
                                   shells_i[jj][0] * a[2]
                                   + shells_i[jj][1] * b[2]
                                   + shells_i[jj][2] * c[2]);
        rvec_i[jj*3 + 0] =
          (*ghost_positions)(j,0) - (*ghost_positions)(i,0) + offset[0];
        rvec_i[jj*3 + 1] =
          (*ghost_positions)(j,1) - (*ghost_positions)(i,1) + offset[1];
        rvec_i[jj*3 + 2] =
          (*ghost_positions)(j,2) - (*ghost_positions)(i,2) + offset[2];
      }
    }
  }
}


double Box::calc_dist(int i, int j, double& dx, double& dy, double& dz) const {
  if (!ghost_positions)
    throw runtime_error("Ghosts no initialized!");
  dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
  dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
  dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
  if (kim_neighbor_mode == KIM_mi_opbc_f) {
    if (abs(dx) > 0.5 * box_side_lengths_[0])
      dx -= (dx / abs(dx)) * box_side_lengths_[0];
    if (abs(dy) > 0.5 * box_side_lengths_[1])
      dy -= (dy / abs(dy)) * box_side_lengths_[1];
    if (abs(dz) > 0.5 * box_side_lengths_[2])
      dz -= (dz / abs(dz)) * box_side_lengths_[2];
  }
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


void Box::deform(Voigt6<double> defmatrix) {
  // Check MI_OPBC_F constraints.
  if (kim_neighbor_mode == KIM_mi_opbc_f
      && (abs(defmatrix(3)) > DELTA
          || abs(defmatrix(4)) > DELTA
          || abs(defmatrix(5)) > DELTA))
    throw runtime_error("cannot shear box with MI_OPBC_F neighbor list.");

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
  update_ghost_rvecs();
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
  switch (kim_neighbor_mode) {
  case KIM_mi_opbc_f:
    // No ghost atoms needed.
    return Vec3D<unsigned>(0, 0, 0);
  case KIM_cluster:
  case KIM_neigh_pure_f:
  case KIM_neigh_rvec_f:
    {
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
    break;
  default:
    throw runtime_error("unknown neighbor list mode.");
  }
}





// Template for the test KIM descriptor file ///////////////////////////

/* Atom types and neighbor list modes will be added dynamically. */
static const string KIM_DESC =
  "Unit_length      := A\n"
  "Unit_energy      := eV\n"
  "Unit_charge      := e\n"
  "Unit_temperature := K\n"
  "Unit_time        := ps\n"
  "\n"
  "ZeroBasedLists    flag\n"
  "\n"
  "# Model input\n"
  "numberOfParticles    integer   none      []\n"
  "numberParticleTypes  integer   none      []\n"
  "particleTypes        integer   none      [numberOfParticles]\n"
  "coordinates          double    length    [numberOfParticles,3]\n"
  "get_neigh            method    none      []\n"
  "neighObject          pointer   none      []\n"
  "boxSideLengths       double    length    [3]\n"
  "\n"
  "# Model output\n"
  "destroy              method    none      []\n"
  "compute              method    none      []\n"
  "reinit               method    none      []\n"
  "cutoff               double    length    []\n"
  "energy               double    energy    []\n"
  "forces               double    force     [numberOfParticles,3]\n"
  "particleEnergy       double    energy    [numberOfParticles]\n"
  "virial               double    energy    [6]\n"
  "particleVirial       double    energy    [numberOfParticles,6]\n"
  ;


// Compute /////////////////////////////////////////////////////////////


Compute::Compute(unique_ptr<Box> box, const string& modelname)
  : box_(move(box)), modelname_(modelname),
    neighbor_mode(box_->kim_neighbor_mode)
{
  int status;

  // Assemble KIM descriptor file.
  string descriptor = KIM_DESC;
  descriptor += "Neigh_LocaAccess  flag\n";
  switch (neighbor_mode) {
  case KIM_cluster:
    descriptor += "CLUSTER           flag\n";
    kim_wants_rvec = false;
    break;
  case KIM_mi_opbc_f:
    descriptor += "MI_OPBC_F         flag\n";
    kim_wants_rvec = false;
    break;
  case KIM_neigh_pure_f:
    descriptor += "NEIGH_PURE_F      flag\n";
    kim_wants_rvec = false;
    break;
  case KIM_neigh_rvec_f:
    descriptor += "NEIGH_RVEC_F      flag\n";
    kim_wants_rvec = true;
    break;
  default:
    throw runtime_error("Unsupported neighbor list mode.");
  }

  // Query and store particle type codes.  For this we first get a
  // dummy model instance that can be queried for all the information
  // in the KIM descriptor.  Also add species to the test descriptor.
  KIM_API_model* query;
  status = KIM_API_model_info(&query, modelname.c_str());
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  char* pt = query->get_model_partcl_typs(&ntypes, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  for (int i = 0; i != ntypes; ++i) {
    const string t(&pt[i*KIM_KEY_STRING_LENGTH]);
    int code = query->get_partcl_type_code(t.c_str(), &status);
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    partcl_type_codes[t] = code;
    partcl_type_names[code] = t;
    // Add to descriptor.
    descriptor += t + " spec " + to_string(code) + "\n";
  }
  free(pt);
  // Destroy the dummy model again.
  KIM_API_free(&query, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Check if the atom types in the box are supported by the model.
  for (int i = 0; i != box_->types.extent(0); ++i)
    if (partcl_type_names.find(box_->types(i)) == partcl_type_names.end())
      // Unknown particle code.
      throw runtime_error("Unknown particle code: "+to_string(box_->types(i)));

  // Init KIM.
  status = KIM_API_string_init(&model, descriptor.c_str(),
                               modelname.c_str());
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Set constant length references to data for KIM.  We only know how
  // long the variable length data is after calculating the ghost
  // atoms.  But we can only do that after initializing the model,
  // which outputs the cutoff.  So we have to make do with these for
  // now.
  model->setm_data(&status, 7*4,
                   "numberOfParticles",   1, &box_->nall,                1,
                   "numberParticleTypes", 1, &ntypes,                    1,
                   "neighObject",         1, this,                       1,
                   "boxSideLengths",      3, &box_->box_side_lengths[0], 1,
                   "cutoff",              1, &cutoff,                    1,
                   "energy",              1, &energy,                    1,
                   "virial",              6, &virial(0),                 1
                   );
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Pass methods to KIM.
  status = model->set_method("get_neigh", 1, (func_ptr) &get_neigh);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Init KIM model.
  status = model->model_init();
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Store free parameter names of the model.
  int n_params;
  char* param_names = model->get_free_params(&n_params, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  for (int i = 0; i != n_params; ++i) {
    const string param_name = &param_names[i*KIM_KEY_STRING_LENGTH];
    const int param_index = model->get_index(param_name.c_str(), &status);
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    const int param_rank = model->get_rank_by_index(param_index, &status);
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    vector<int> param_shape_(param_rank);
    model->get_shape_by_index(param_index, &param_shape_[0], &status);
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    vector<unsigned> param_shape;
    for (int j = 0; j != param_rank; ++j)
      param_shape.push_back(param_shape_[j]);
    FreeParam fp(param_name, param_shape, param_index);
    free_parameters.push_back(fp);
    // Fancy way to do free_parameter_map[param_name] = fp; but
    // without assignment or copy.
    free_parameter_map.emplace(piecewise_construct,
                               forward_as_tuple(param_name),
                               forward_as_tuple(param_name, param_shape,
                                                param_index));
  }
  free(param_names);

  // Now that we know the cutoff, we calculate neighbor lists and
  // ghost atoms.  Then we allocate memory for variable length data
  // and pass it to KIM.
  update_neighbor_list();
}

Compute::~Compute() {
  int status;
  status = model->model_destroy();
  if (status < KIM_STATUS_OK)
    cout << string("KIM error in line ") + to_string(__LINE__)
      + string(" of file ") + string(__FILE__) << endl;
  KIM_API_free(&model, &status);
  if (status < KIM_STATUS_OK)
    cout << string("KIM error in line ") + to_string(__LINE__)
      + string(" of file ") + string(__FILE__) << endl;
}


void Compute::compute() {
  const int status = model->model_compute();
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

  // Tricks to fix values when ghost atoms are in the game.
  if (box_->nghosts)
    switch (neighbor_mode) {
    case KIM_cluster:
      // Need to calculate whole box values from particle values due to
      // the fact that ghost atoms were included in the summation.
      energy = 0.0;
      virial(0) = 0.0; virial(1) = 0.0; virial(2) = 0.0;
      virial(3) = 0.0; virial(4) = 0.0; virial(5) = 0.0;
      for (unsigned i = 0; i != box_->natoms; ++i) {
        energy += particleEnergy[i];
        virial(0) += particleVirial[i*6 + 0];
        virial(1) += particleVirial[i*6 + 1];
        virial(2) += particleVirial[i*6 + 2];
        virial(3) += particleVirial[i*6 + 3];
        virial(4) += particleVirial[i*6 + 4];
        virial(5) += particleVirial[i*6 + 5];
      }
      break;
    case KIM_mi_opbc_f:
    case KIM_neigh_rvec_f:
      // These have no ghost atoms.
      throw runtime_error("this should never happen: there are ghost atoms "
                          "although we use MI_OPBC_F.");
    case KIM_neigh_pure_f:
      for (unsigned i = box_->natoms; i != box_->nall; ++i) {
        const unsigned central = i % box_->natoms;
        particleEnergy[central] += particleEnergy[i];
        forces[3*central + 0] += forces[3*i + 0];
        forces[3*central + 1] += forces[3*i + 1];
        forces[3*central + 2] += forces[3*i + 2];
        particleVirial[6*central + 0] += particleVirial[6*i + 0];
        particleVirial[6*central + 1] += particleVirial[6*i + 1];
        particleVirial[6*central + 2] += particleVirial[6*i + 2];
        particleVirial[6*central + 3] += particleVirial[6*i + 3];
        particleVirial[6*central + 4] += particleVirial[6*i + 4];
        particleVirial[6*central + 5] += particleVirial[6*i + 5];
      }
      break;
    default:
      throw runtime_error("unsupported neighbor mode.");
  }

  /*
  cout << "                                                                 "
       << "               "
       << forces[0] << " " << forces[1] << " " << forces[2]
       << "  ::  "
       << particleVirial[0] << " " << particleVirial[1] << " "
       << particleVirial[2] << " " << particleVirial[3] << " "
       << particleVirial[4] << " " << particleVirial[5] << " "
       << "  ::  "
       << particleEnergy[0]
       << endl;
  */
}


int Compute::get_neigh(KIM_API_model** kimmdl,
                       const int *mode, const int* request,
                       int *particle, int *numnei, int **nei1particle,
                       double **rij) {
  KIM_API_model& model = **kimmdl;
  int status;
  Compute& c = *( (Compute*)model.get_data("neighObject", &status) );
  if (status < KIM_STATUS_OK)
    return KIM_STATUS_FAIL;

  // Get central atom.
  int i;
  if (*mode == 0) { // Iterator mode.
    if (*request == 0) { // Reset requested.
      c.kim_iter_pos = 0;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    } else {             // Increment requested.
      i = c.kim_iter_pos;
      // Cancel if all non-ghost atoms exhausted.
      if (i >= static_cast<int>(c.box_->natoms))
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      ++c.kim_iter_pos;
    }
  } else {          // Locator mode.
    i = *request;
  }

  // Prepare return values.
  if (i < static_cast<int>(c.box_->natoms)) {
    // Not a ghost atom.
    const vector<int>& neighbors = c.box_->get_neighbors(i);
    *particle = i;
    *numnei = neighbors.size();
    // We can get away with the const casts because (at least by
    // convention) KIM guarantees that the model will not fiddle with
    // these arrays.
    *nei1particle = const_cast<int*>(&neighbors[0]);
    *rij = c.kim_wants_rvec
      ? const_cast<double*>(c.box_->get_neighbor_rvecs_ptr(i))
      : nullptr;
  } else {
    // Ghost atom.
    *particle = i;
    *numnei = 0;
    *nei1particle = nullptr;
    *rij = nullptr;
  }

  return KIM_STATUS_OK;
}

void Compute::update_neighbor_list() {
  int status;
  // TODO: For now the skin is 0.
  const bool arrays_changed = box_->update_neighbor_list(cutoff, 0.0);
  if (arrays_changed) {
    const unsigned n = box_->nall;
    forces.resize(3 * n);
    particleEnergy.resize(n);
    particleVirial.resize(6 * n);
    model->setm_data(&status, 5*4, // TODO: use indices!
                     "coordinates",    3*n, box_->get_positions_ptr(), 1,
                     "particleTypes",    n, box_->get_types_ptr(),     1,
                     "forces",         3*n, &forces[0],                1,
                     "particleEnergy",   n, &particleEnergy[0],        1,
                     "particleVirial", 6*n, &particleVirial[0],        1
                     );
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
  }
}

template<typename T>
void Compute::set_parameter_impl(const string& param_name,
                                 const vector<unsigned>& indices,
                                 T value, bool reinit) {
  // Get parameter data.
  const auto it = free_parameter_map.find(param_name);
  if (it == free_parameter_map.end())
    throw runtime_error("Unknown free parameter: " + param_name);
  const FreeParam& param = it->second;
  // Check rank.
  if (indices.size() != param.rank)
    throw runtime_error("Rank of indices ("
                        + to_string(indices.size())
                        + ") doesn't match rank of parameter ("
                        + to_string(param.rank) + ")");
  // Set it.
  int status;
  T* p = (T*)model->get_data_by_index(param.kim_index, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  unsigned index = 0;
  for (unsigned k = 0; k != param.rank; ++k) {
    unsigned product = 1;
    for (unsigned l = k + 1; l != param.rank; ++l)
      product *= param.shape[l];
    index += product * indices[k];
  }
  p[index] = value;

  // Reinit model.
  if (reinit) {
    status = model->model_reinit();
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
  }
}

void Compute::set_parameter(const string& param_name,
                            const vector<unsigned>& indices,
                            double value, bool reinit) {
  // TODO: check if this uses the correct type.     
  set_parameter_impl<double>(param_name, indices, value, reinit);
}

void Compute::set_parameter(const string& param_name,
                            const vector<unsigned>& indices,
                            int value, bool reinit) {
  // TODO: check if this uses the correct type.     
  set_parameter_impl<int>(param_name, indices, value, reinit);
}


// Fitting /////////////////////////

double Compute::obj_func_box(const vector<double>& x, vector<double>& grad,
                             void* f_data) {
  Compute& c = *(static_cast<Compute*>(f_data));
  if (x.size() != 3)
    throw runtime_error("This objective function only works with "
                        "3 parameters.");
  // Scale box, we assume the neighbor list doesn't change (TODO?).
  c.box_->scale_to(x[0], x[1], x[2]);
  c.compute();
  ++c.fit_counter;
  // Gradient.
  if (!grad.empty())
    throw runtime_error("Gradient not supported for box optimization.");
  // Objective function value is energy.
  return c.energy;
}

double Compute::optimize_box(double ftol_abs, unsigned maxeval) {
  vector<double> lengths = { box_->box_side_lengths[0],
                             box_->box_side_lengths[1],
                             box_->box_side_lengths[2] };
  vector<double> lb = { 0.0, 0.0, 0.0 }; // Lengths may not be negative.
  // TODO: We use a gradient-free algorithm for now, although the
  // gradient should be obtainable from the virial (is it even exactly
  // the virial?).
  double obj_val;
  nlopt::opt optimizer(nlopt::LN_SBPLX, 3);
  optimizer.set_min_objective(Compute::obj_func_box, this);
  optimizer.set_lower_bounds(lb);
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
  c.box_->update_ghost_rvecs();
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
  // Get coordinates as vector.
  Array2D<double>& pos = box_->positions;
  const int nparams = pos.extent(0)*pos.extent(1);
  const double* p = &pos(0,0);
  vector<double> pos_vec(p, p+nparams);
  // Other good algorithms:
  //   * NLOPT_LD_LBFGS
  //   * NLOPT_LD_TNEWTON_PRECOND_RESTART
  //   * NLOPT_LD_VAR2
  double obj_val;
  nlopt::opt optimizer(nlopt::LD_TNEWTON_PRECOND_RESTART, nparams);
  optimizer.set_min_objective(Compute::obj_func_pos, this);
  fit_counter = 0;
  optimizer.set_maxeval(maxeval);
  optimizer.set_ftol_abs(ftol_abs);
  optimizer.optimize(pos_vec, obj_val);
  return obj_val;
}



void Compute::switch_boxes(Compute& other) {
  // Check neighbor list mode.
  if (neighbor_mode != other.neighbor_mode)
    throw runtime_error("Neighbor list modes do not match.");
  // Check model.
  if (modelname_ != other.modelname_)
    throw runtime_error("KIM models do not match.");
  // Swap boxes.
  box_.swap(other.box_);
  // Reset KIM data on this object.
  int status;
  const unsigned n = box_->nall;
  forces.resize(3 * n);
  particleEnergy.resize(n);
  particleVirial.resize(6 * n);
  model->setm_data(&status, 7*4, // TODO: use indices!
                   "numberOfParticles",   1, &box_->nall,                1,
                   "boxSideLengths",      3, &box_->box_side_lengths[0], 1,
                   "coordinates",       3*n, box_->get_positions_ptr(),  1,
                   "particleTypes",       n, box_->get_types_ptr(),      1,
                   "forces",            3*n, &forces[0],                 1,
                   "particleEnergy",      n, &particleEnergy[0],         1,
                   "particleVirial",    6*n, &particleVirial[0],         1
                   );
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  // Reset KIM data on other object.
  const unsigned n2 = other.box_->nall;
  Box& b = *other.box_;
  other.forces.resize(3 * n2);
  other.particleEnergy.resize(n2);
  other.particleVirial.resize(6 * n2);
  other.model->setm_data(&status, 7*4, // TODO: use indices!
                         "numberOfParticles", 1, &b.nall,                  1,
                         "boxSideLengths",    3, &b.box_side_lengths[0],   1,
                         "coordinates",    3*n2, b.get_positions_ptr(),    1,
                         "particleTypes",    n2, b.get_types_ptr(),        1,
                         "forces",         3*n2, &other.forces[0],         1,
                         "particleEnergy",   n2, &other.particleEnergy[0], 1,
                         "particleVirial", 6*n2, &other.particleVirial[0], 1
                         );
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
}







void do_something(const char* lat, double lat_const, KIMNeigh neighmode) {
  string neighmode_str;
  switch (neighmode) {
  case KIM_cluster:
    neighmode_str = "CLUSTER";
    break;
  case KIM_mi_opbc_f:
    neighmode_str = "MI_OPBC_F";
    break;
  case KIM_neigh_pure_f:
    neighmode_str = "NEIGH_PURE_F";
    break;
  case KIM_neigh_rvec_f:
    neighmode_str = "NEIGH_RVEC_F";
    break;
  default:
    throw runtime_error("unknown neighbor mode.");
  }
  // Init.
  const vector<int> types{ 0 };
  Compute comp(make_unique<Box>(lat, lat_const*0.9, true, 3, 3, 3,
                                true, true, true,
                                types, neighmode, "box"),
               "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");

  comp.set_parameter("PARAM_FREE_Rc", vector<unsigned>{ 0, 0, 0 }, 2.9, false);
  comp.set_parameter("PARAM_FREE_Dc", vector<unsigned>{ 0, 0, 0 }, 0.15);
  comp.update_neighbor_list();
  /*
  comp.compute();
  cout << comp.get_energy_per_atom() << " eV/atom";
  for (unsigned i = 0; i != 6; ++i)
    printf("  %8.4g", comp.get_virial()(i));
  cout << endl;
  */

  chrono::steady_clock::time_point start;
  chrono::steady_clock::time_point stop;
  chrono::duration<double> delta;

  comp.move_atom(0, 0.2, 0.0, 0.0);

  /*
  cout << comp.get_particle_type_code("Si") << " : "
       << comp.get_particle_type_name(0) << endl;
  */
  printf("== Box scaling ================ %12s ============\n",
         neighmode_str.c_str());
  comp.compute();
  cout << "\nBefore optimization: "
       << comp.get_energy_per_atom() << " eV/atom" << endl;
  try {
    start = chrono::steady_clock::now();
    comp.optimize_box(0.001, 10000);
    stop = chrono::steady_clock::now();
    delta = chrono::duration_cast<chrono::duration<double>>(stop - start);
    cout << "After " << comp.get_optimization_steps() << " steps: "
         << comp.get_energy_per_atom() << " eV/atom"
         << "                       took " << delta.count() << " s."
         << endl;
    cout << "a = " << comp.get_a().abs()/3 << "  "
         << "b = " << comp.get_b().abs()/3 << "  "
         << "c = " << comp.get_c().abs()/3 << "  "
         << "goal = " << lat_const << endl;
  } catch (const exception& e) {
    cout << "NLopt failed :-(" << endl;
  }
  cout << endl;

  printf("== Position scaling =========== %12s ============\n",
         neighmode_str.c_str());
  comp.compute();
  cout << "\nBefore optimization: "
       << comp.get_energy_per_atom() << " eV/atom" << endl;
  try {
    start = chrono::steady_clock::now();
    comp.optimize_positions(0.001, 10000);
    stop = chrono::steady_clock::now();
    delta = chrono::duration_cast<chrono::duration<double>>(stop - start);
    cout << "After " << comp.get_optimization_steps() << " steps: "
         << comp.get_energy_per_atom() << " eV/atom"
         << "                       took " << delta.count() << " s."
         << endl;
  } catch (const exception& e) {
    cout << "NLopt failed :-(" << endl;
  }
  cout << endl;

  printf("== Box scaling II ============= %12s ============\n",
         neighmode_str.c_str());
  comp.compute();
  cout << "\nBefore optimization: "
       << comp.get_energy_per_atom() << " eV/atom" << endl;
  try {
    start = chrono::steady_clock::now();
    comp.optimize_box(0.001, 10000);
    stop = chrono::steady_clock::now();
    delta = chrono::duration_cast<chrono::duration<double>>(stop - start);
    cout << "After " << comp.get_optimization_steps() << " steps: "
         << comp.get_energy_per_atom() << " eV/atom"
         << "                       took " << delta.count() << " s."
         << endl;
    cout << "a = " << comp.get_a().abs()/3 << "  "
         << "b = " << comp.get_b().abs()/3 << "  "
         << "c = " << comp.get_c().abs()/3 << "  "
         << "goal = " << lat_const << endl;
  } catch (const exception& e) {
    cout << "NLopt failed :-(" << endl;
  }
  cout << endl;
}


/*
int main() {
  const char lat[] = "fcc";
  const double lat_const = 3.940;
  vector<int> types;
  types.push_back(0);
  //types.push_back(0);
  Box b(lat, lat_const, true, 2, 2, 2, true, true, true,
        types, KIM_neigh_pure_f, "box");

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
  cout << endl << endl;

  // b.write_to("dump");

  //////////////////////////////////////////////////////////////////////
  do_something(lat, lat_const, KIM_mi_opbc_f);
  do_something(lat, lat_const, KIM_neigh_pure_f);
  do_something(lat, lat_const, KIM_neigh_rvec_f);
  return 0;
  do_something(lat, lat_const, KIM_cluster);
  do_something(lat, lat_const, KIM_mi_opbc_f);
  do_something(lat, lat_const, KIM_neigh_pure_f);
  do_something(lat, lat_const, KIM_neigh_rvec_f);

  return 0;
}
*/
