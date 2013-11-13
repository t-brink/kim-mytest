#include <iostream>
#include <cstdio>
#include <utility>
#include <cmath>

#include <KIM_API_C.h>
#include <KIM_API_status.h>

#include "mytest.hpp"

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
    neigh_rvec_(neighmode == KIM_neigh_rvec_f ? natoms_ : 0)
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
  if (neighmode == KIM_neigh_rvec_f)
    neigh_rvec_.resize(natoms_);
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
            if (kim_neighbor_mode == KIM_neigh_rvec_f) {
              neigh_rvec_[i].push_back(dx);
              neigh_rvec_[i].push_back(dy);
              neigh_rvec_[i].push_back(dz);
              neigh_rvec_[j].push_back(-dx);
              neigh_rvec_[j].push_back(-dy);
              neigh_rvec_[j].push_back(-dz);
            }
          }
        }
        // Iterate over ghosts.
        for (unsigned j = natoms_; j != nall_; ++j) {
          const double dx = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
          const double dy = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
          const double dz = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
          if (dx*dx + dy*dy + dz*dz < cutsq) {
            neigh_list_[i].push_back(j);
            if (kim_neighbor_mode == KIM_neigh_rvec_f) {
              neigh_rvec_[i].push_back(dx);
              neigh_rvec_[i].push_back(dy);
              neigh_rvec_[i].push_back(dz);
            }
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
            if (kim_neighbor_mode == KIM_neigh_rvec_f) {
              neigh_rvec_[i].push_back(dx);
              neigh_rvec_[i].push_back(dy);
              neigh_rvec_[i].push_back(dz);
              neigh_rvec_[j].push_back(-dx);
              neigh_rvec_[j].push_back(-dy);
              neigh_rvec_[j].push_back(-dz);
            }
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


void Box::update_ghost_rvecs() {
  update_ghosts();
  if (kim_neighbor_mode == KIM_neigh_rvec_f)
    for (unsigned i = 0; i != natoms_; ++i) {
      const unsigned jjmax = neigh_list_[i].size();
      const auto& neigh_i = neigh_list_[i];
      auto& rvec_i = neigh_rvec_[i];
      for (unsigned jj = 0; jj != jjmax; ++jj) {
        const unsigned j = neigh_i[jj];
        rvec_i[jj*3 + 0] = (*ghost_positions)(j,0) - (*ghost_positions)(i,0);
        rvec_i[jj*3 + 1] = (*ghost_positions)(j,1) - (*ghost_positions)(i,1);
        rvec_i[jj*3 + 2] = (*ghost_positions)(j,2) - (*ghost_positions)(i,2);
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
  : box_(move(box)), neighbor_mode(box_->kim_neighbor_mode)
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
    // Add to descriptor.
    descriptor += t + " spec " + to_string(code) + "\n";
  }
  // Destroy the dummy model again.
  KIM_API_free(&query, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));

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

  // Now that we know the cutoff, we calculate neighbor lists and
  // ghost atoms.  TODO: For now the skin is 0.
  box_->update_neighbor_list(cutoff, 0.0);

  // Allocate memory for variable length data and pass it to KIM.
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








int main() {
  vector<int> types;
  types.push_back(0);
  //types.push_back(0);
  Box b("sc", 2.525, true, 3, 3, 3, true, true, true,
        types, KIM_mi_opbc_f, "box");

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

  // b.write_to("dump");

  //////////////////////////////////////////////////////////////////////

  {
    Compute comp(make_unique<Box>("sc", 2.525, true, 3, 3, 3,
                                  true, true, true,
                                  types, KIM_cluster, "box"),
                 "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
    comp.compute();
    cout << comp.get_energy_per_atom() << " eV/atom";
    for (unsigned i = 0; i != 6; ++i)
      printf("  %8.4g", comp.get_virial()(i));
    cout << endl;
  }

  {
    Compute comp(make_unique<Box>("sc", 2.525, true, 3, 3, 3,
                                  true, true, true,
                                  types, KIM_mi_opbc_f, "box"),
                 "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
    comp.compute();
    cout << comp.get_energy_per_atom() << " eV/atom";
    for (unsigned i = 0; i != 6; ++i)
      printf("  %8.4g", comp.get_virial()(i));
    cout << endl;
  }

  {
    Compute comp(make_unique<Box>("sc", 2.525, true, 3, 3, 3,
                                  true, true, true,
                                  types, KIM_neigh_pure_f, "box"),
                 "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
    comp.compute();
    cout << comp.get_energy_per_atom() << " eV/atom";
    for (unsigned i = 0; i != 6; ++i)
      printf("  %8.4g", comp.get_virial()(i));
    cout << endl;
  }

  {
    Compute comp(make_unique<Box>("sc", 2.525, true, 3, 3, 3,
                                  true, true, true,
                                  types, KIM_neigh_rvec_f, "box"),
                 "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000");
    comp.compute();
    cout << comp.get_energy_per_atom() << " eV/atom";
    for (unsigned i = 0; i != 6; ++i)
      printf("  %8.4g", comp.get_virial()(i));
    cout << endl;
  }

  return 0;
}
