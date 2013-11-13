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

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <memory>
#include <map>
//#include <algorithm>

#include <KIM_API.h>
#include <KIM_API_C.h>
#include <KIM_API_status.h>

#include <nlopt.hpp>

#include "ndarray.hpp"

using namespace std;
using namespace mytest;


// Constants ///////////////////////////////////////////////////////////

const char MODELNAME[] = "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_000";


// Template for the test KIM descriptor file ///////////////////////////

/* Atom types and neighbor list modes will be added dynamically. */
const string KIM_DESC =
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


// Utils ///////////////////////////////////////////////////////////////
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique(Args&& ...args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/// Modulo with the same semantics as python (or maths).
double pmod(double x, double N) {
  return (x < 0) ? fmod(fmod(x, N) + N, N) : fmod(x, N);
}

// Data structures /////////////////////////////////////////////////////


enum KIMNeigh {
  KIM_cluster, KIM_mi_opbc_f, KIM_neigh_pure_f, KIM_neigh_rvec_f
};

class Box {
protected:
  Box(int n, double lenx, double leny, double lenz, const string& lattice_type)
    : N(n), // Number of atoms.
      nghosts(0), // No ghosts, yet
      box_size_(lenx, leny, lenz),
      coordinates_(N,3),
      types_(N),
      neigh_list(N),
      neigh_rvec(N),
      lattice_type_(lattice_type)
  {}

public:
  virtual ~Box() {}

  /*! Return number of additional boxes in a, b, and c direction. */
  Vec3D<unsigned> nshells(double cutoff) const { // TODO: PBC on/off  
    // TODO: once we have general vectors use them.  
    const Vec3D<double> a(box_size_[0], 0.0, 0.0);
    const Vec3D<double> b(0.0, box_size_[1], 0.0);
    const Vec3D<double> c(0.0, 0.0, box_size_[2]);
    const Vec3D<double> bc = cross(b, c);
    const Vec3D<double> ac = cross(a, c);
    const Vec3D<double> ab = cross(a, b);
    const double V = abs(dot(a, bc));
    const double da = V / bc.abs();
    const double db = V / ac.abs();
    const double dc = V / ab.abs();
    return Vec3D<unsigned>(ceil(cutoff/da),
                           ceil(cutoff/db),
                           ceil(cutoff/dc));
  }

  /*! Update ghost positions from central cell.

    Assumes that the ghost coordinate array has the right size.
  */
  void update_ghosts(const Vec3D<unsigned>& shells) {
    // TODO: once we have general vectors use them.  
    const Vec3D<double> a(box_size_[0], 0.0, 0.0);
    const Vec3D<double> b(0.0, box_size_[1], 0.0);
    const Vec3D<double> c(0.0, 0.0, box_size_[2]);
    unsigned ii = 0;
    for (unsigned i = 0; i != N; ++i) {
      (*ghost_coordinates_)(i, 0) = coordinates_(i, 0);
      (*ghost_coordinates_)(i, 1) = coordinates_(i, 1);
      (*ghost_coordinates_)(i, 2) = coordinates_(i, 2);
      ++ii;
    }
    const int alo = -static_cast<int>(shells[0]);
    const int ahi = shells[0];
    const int blo = -static_cast<int>(shells[1]);
    const int bhi = shells[1];
    const int clo = -static_cast<int>(shells[2]);
    const int chi = shells[2];
    for (int aa = alo; aa <= ahi; ++aa)
    for (int bb = blo; bb <= bhi; ++bb)
    for (int cc = clo; cc <= chi; ++cc) {
      if (aa == 0 && bb == 0 && cc == 0) continue;
      for (unsigned i = 0; i != N; ++i) {
        const Vec3D<double> offset = 1.0*aa*a + 1.0*bb*b + 1.0*cc*c;
        (*ghost_coordinates_)(ii, 0) = coordinates_(i, 0) + offset[0];
        (*ghost_coordinates_)(ii, 1) = coordinates_(i, 1) + offset[1];
        (*ghost_coordinates_)(ii, 2) = coordinates_(i, 2) + offset[2];
        ghost_types_[ii] = types_[i];
        ++ii;
      }
    }
  }

  void init_neighbor_list(double cutoff) {
    const Vec3D<unsigned> shells = nshells(cutoff);
    const unsigned new_nghosts =
      N * ((2*shells[0]+1) * (2*shells[1]+1) * (2*shells[2]+1) - 1);
    if (new_nghosts != nghosts) {
      nghosts = new_nghosts;
      ghost_coordinates_ = make_unique< Array2D<double> >(N+nghosts, 3);
      ghost_types_.resize(N+nghosts);
      update_ghosts(shells);
    }
    // Clear neighbor lists.
    for (auto& l : neigh_list)
      l.clear();
    for (auto& r : neigh_rvec)
      r.clear();
    // Fill neighbor lists again.
    const double cutsq = cutoff * cutoff;
    for (unsigned i = 0; i != N; ++i) {
      // Go over neighbors in the central cell.
      for (unsigned j = i+1; j != N; ++j) {
        const double dx = coordinates_(j,0) - coordinates_(i,0);
        const double dy = coordinates_(j,1) - coordinates_(i,1);
        const double dz = coordinates_(j,2) - coordinates_(i,2);
        if (dx*dx + dy*dy + dz*dz < cutsq) {
          neigh_list[i].push_back(j);
          neigh_list[j].push_back(i);
          neigh_rvec[i].push_back(dx);
          neigh_rvec[i].push_back(dy);
          neigh_rvec[i].push_back(dz);
          neigh_rvec[j].push_back(-dx);
          neigh_rvec[j].push_back(-dy);
          neigh_rvec[j].push_back(-dz);
        }
      }
      // Iterate over ghosts.
      for (unsigned j = N; j != N+nghosts; ++j) {
        const double dx = (*ghost_coordinates_)(j,0) - coordinates_(i,0);
        const double dy = (*ghost_coordinates_)(j,1) - coordinates_(i,1);
        const double dz = (*ghost_coordinates_)(j,2) - coordinates_(i,2);
        if (dx*dx + dy*dy + dz*dz < cutsq) {
          neigh_list[i].push_back(j);
          neigh_rvec[i].push_back(dx);
          neigh_rvec[i].push_back(dy);
          neigh_rvec[i].push_back(dz);
        }
      }
    }
  }

  const unsigned& natoms() const {
    return N;
  }

  unsigned natoms_all() const {
    return N+nghosts;
  }

  const Vec3D<double>& box_size() const {
    return box_size_;
  }

  double get_volume() const {
    return box_size_[0] * box_size_[1] * box_size_[2];
  }

  Array2D<double>& coordinates() {
    return coordinates_;
  }

  const Array2D<double>& ghost_coordinates() const {
    return *ghost_coordinates_;
  }

  const vector<int>& types() const {
    return types_;
  }

  const vector<int>& ghost_types() const {
    return ghost_types_;
  }

  const string& lattice_type() const {
    return lattice_type_;
  }

  // Distance using PBC (minimum image convention).
  double dist(int i, int j) const {
    double dx, dy, dz;
    return dist(i, j, dx, dy, dz);
  }

  double dist(int i, int j, double& dx, double& dy, double& dz) const {
    dx = (*ghost_coordinates_)(j,0) - coordinates_(i,0);
    dy = (*ghost_coordinates_)(j,1) - coordinates_(i,1);
    dz = (*ghost_coordinates_)(j,2) - coordinates_(i,2);
    /*
    if (abs(dx) > 0.5 * box_size_[0])
      dx -= (dx/abs(dx))*box_size_[0];
    if (abs(dy) > 0.5 * box_size_[1])
      dy -= (dy/abs(dy))*box_size_[1];
    if (abs(dz) > 0.5 * box_size_[2])
      dz -= (dz/abs(dz))*box_size_[2];
    */
    return sqrt(dx*dx + dy*dy + dz*dz);
  }

  const vector< vector<int> >& get_neigh_list() const {
    return neigh_list;
  }

  const vector< vector<double> >& get_neigh_rvec() const {
    return neigh_rvec;
  }

  /// Write LAMMPS dump file.
  const void write_to(const string& filename) const {
    ofstream myfile;
    myfile.open(filename);
    myfile << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"
           << N
           << "\nITEM: BOX BOUNDS pp pp pp\n"
           << 0 << " " << box_size_[0] << "\n"
           << 0 << " " << box_size_[1] << "\n"
           << 0 << " " << box_size_[2] << "\n"
           << "ITEM: ATOMS id type x y z\n";
    for (unsigned i = 0; i != N; ++i)
      myfile << i << " " << types_[i] << " "
             << coordinates_(i,0) << " "
             << coordinates_(i,1) << " "
             << coordinates_(i,2) << "\n";
  }

  const void write_to2(const string& filename) const {
    ofstream myfile;
    myfile.open(filename);
    myfile << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"
           << N+nghosts
           << "\nITEM: BOX BOUNDS pp pp pp\n"
           << 0 << " " << box_size_[0] << "\n"
           << 0 << " " << box_size_[1] << "\n"
           << 0 << " " << box_size_[2] << "\n"
           << "ITEM: ATOMS id type x y z\n";
    for (unsigned i = 0; i != N; ++i)
      myfile << i << " " << types_[i] << " "
             << coordinates_(i,0) << " "
             << coordinates_(i,1) << " "
             << coordinates_(i,2) << "\n";
    for (unsigned i = N; i != N+nghosts; ++i)
      myfile << i << " " << 2 << " "
             << (*ghost_coordinates_)(i,0) << " "
             << (*ghost_coordinates_)(i,1) << " "
             << (*ghost_coordinates_)(i,2) << "\n";
  }

  /// Scale box by a factor.
  void scale(double factor_x, double factor_y, double factor_z) {
    box_size_[0] *= factor_x;
    box_size_[1] *= factor_y;
    box_size_[2] *= factor_z;
    for (unsigned i = 0; i != N; ++i) {
      coordinates_(i,0) *= factor_x;
      coordinates_(i,1) *= factor_y;
      coordinates_(i,2) *= factor_z;
    }
  }

  void scale(double factor) {
    box_size_ *= factor;
    coordinates_ *= factor;
  };

  /// Scale box so it conforms to given size.
  void scale_to(double x, double y, double z) {
    const double factor_x = x / box_size_[0];
    const double factor_y = y / box_size_[1];
    const double factor_z = z / box_size_[2];
    box_size_[0] = x;
    box_size_[1] = y;
    box_size_[2] = z;
    for (unsigned i = 0; i != N; ++i) {
      coordinates_(i,0) *= factor_x;
      coordinates_(i,1) *= factor_y;
      coordinates_(i,2) *= factor_z;
    }
  }

protected:
  const unsigned N; // Number of atoms.
  unsigned nghosts; // Number of ghost atoms.
  Vec3D<double> box_size_;
  Array2D<double> coordinates_;
  unique_ptr< Array2D<double> > ghost_coordinates_;
  vector<int> types_;
  vector<int> ghost_types_;
  vector< vector<int> > neigh_list;
  vector< vector<double> > neigh_rvec; // The inner vector is actually
                                       // an N*3 array.
  const string lattice_type_;
};


class dimerBox : public Box {
public:
  dimerBox(double r0)
    : Box(2, 100.0, 100.0, 100.0, "dimer") // 2 atoms -> dimer
  {
    coordinates_(0, 0) = 50.0 - r0/2;
    coordinates_(0, 1) = 50.0;
    coordinates_(0, 2) = 50.0;
    coordinates_(1, 0) = 50.0 + r0/2;
    coordinates_(1, 1) = 50.0;
    coordinates_(1, 2) = 50.0;
    for (unsigned i = 0; i != N; ++i)
      types_[i] = 0;//i % 2; // Si and C mixed    TODO      
  }
};


class scBox : public Box {
public:
  scBox(double a, int rx, int ry, int rz)
    : Box(rx*ry*rz, a*rx, a*ry, a*rz, "sc") // 1 atom in sc unit cell
  {
    unsigned i = 0;
    for (int x = 0; x != rx; ++x)
      for (int y = 0; y != ry; ++y)
        for (int z = 0; z != rz; ++z) {
          coordinates_(i, 0) = x*a; // 0, 0, 0
          coordinates_(i, 1) = y*a;
          coordinates_(i, 2) = z*a;
          i += 1;
        }
    for (i = 0; i != N; ++i)
      types_[i] = 0;//i % 2; // Si and C mixed    TODO      
  }
};


class bccBox : public Box {
public:
  bccBox(double a, int rx, int ry, int rz)
    : Box(2*rx*ry*rz, a*rx, a*ry, a*rz, "bcc") // 2 atoms in fcc unit cell
  {
    const double a2 = a / 2.0;
    unsigned i = 0;
    for (int x = 0; x != rx; ++x)
      for (int y = 0; y != ry; ++y)
        for (int z = 0; z != rz; ++z) {
          coordinates_(i, 0) = x*a; // 0, 0, 0
          coordinates_(i, 1) = y*a;
          coordinates_(i, 2) = z*a;
          coordinates_(i+1, 0) = a2 + x*a; // a/2, a/2, a/2
          coordinates_(i+1, 1) = a2 + y*a;
          coordinates_(i+1, 2) = a2 + z*a;
          i += 2;
        }
    for (i = 0; i != N; ++i)
      types_[i] = 0;//i % 2; // Si and C mixed    TODO      
  }
};


class fccBox : public Box {
public:
  fccBox(double a, int rx, int ry, int rz)
    : Box(4*rx*ry*rz, a*rx, a*ry, a*rz, "fcc") // 4 atoms in fcc unit cell
  {
    const double a2 = a / 2.0;
    unsigned i = 0;
    for (int x = 0; x != rx; ++x)
      for (int y = 0; y != ry; ++y)
        for (int z = 0; z != rz; ++z) {
          coordinates_(i, 0) = x*a; // 0, 0, 0
          coordinates_(i, 1) = y*a;
          coordinates_(i, 2) = z*a;
          coordinates_(i+1, 0) =  0 + x*a; // 0, a/2, a/2
          coordinates_(i+1, 1) = a2 + y*a;
          coordinates_(i+1, 2) = a2 + z*a;
          coordinates_(i+2, 0) = a2 + x*a; // a/2, 0, a/2
          coordinates_(i+2, 1) =  0 + y*a;
          coordinates_(i+2, 2) = a2 + z*a;
          coordinates_(i+3, 0) = a2 + x*a; // a/2, a/2, 0
          coordinates_(i+3, 1) = a2 + y*a;
          coordinates_(i+3, 2) =  0 + z*a;
          i += 4;
        }
    for (i = 0; i != N; ++i)
      types_[i] = 0;//i % 2; // Si and C mixed    TODO      
  }
};


class diamondBox : public Box {
public:
  diamondBox(double a, int rx, int ry, int rz)
    : Box(8*rx*ry*rz, a*rx, a*ry, a*rz, "diamond") // 8 atoms in unit cell
  {
    const double a2 = a / 2.0;
    const double a4 = a / 4.0;
    const double a34 = 0.75 * a;
    unsigned i = 0;
    for (int x = 0; x != rx; ++x)
      for (int y = 0; y != ry; ++y)
        for (int z = 0; z != rz; ++z) {
          coordinates_(i, 0) = x*a; // 0, 0, 0
          coordinates_(i, 1) = y*a;
          coordinates_(i, 2) = z*a;
          coordinates_(i+1, 0) =  0  + x*a; // 0, a/2, a/2
          coordinates_(i+1, 1) = a2  + y*a;
          coordinates_(i+1, 2) = a2  + z*a;
          coordinates_(i+2, 0) = a2  + x*a; // a/2, 0, a/2
          coordinates_(i+2, 1) =  0  + y*a;
          coordinates_(i+2, 2) = a2  + z*a;
          coordinates_(i+3, 0) = a2  + x*a; // a/2, a/2, 0
          coordinates_(i+3, 1) = a2  + y*a;
          coordinates_(i+3, 2) =  0  + z*a;
          coordinates_(i+4, 0) = a4  + x*a; // a/4, a/4, a/4
          coordinates_(i+4, 1) = a4  + y*a;
          coordinates_(i+4, 2) = a4  + z*a;
          coordinates_(i+5, 0) = a4  + x*a; // a/4, a*3/4, a*3/4
          coordinates_(i+5, 1) = a34 + y*a;
          coordinates_(i+5, 2) = a34 + z*a;
          coordinates_(i+6, 0) = a34 + x*a; // a*3/4, a/4, a*3/4
          coordinates_(i+6, 1) = a4  + y*a;
          coordinates_(i+6, 2) = a34 + z*a;
          coordinates_(i+7, 0) = a34 + x*a; // a*3/4, a*3/4, a/4
          coordinates_(i+7, 1) = a34 + y*a;
          coordinates_(i+7, 2) = a4  + z*a;
          i += 8;
        }
    for (i = 0; i != N; ++i)
      types_[i] = 0;//i % 2; // Si and C mixed    TODO      
  }
};


class NeighborListIterator {
public:
  NeighborListIterator(const Box& box, bool rvec)
    : i(0), b(box), has_rvec(rvec)
  {}
  unsigned i;
  const Box& b;
  bool has_rvec;
};


// Create a box. ///////////////////////////////////////////////////////

unique_ptr<Box> make_box(const string& lattice, double a,
                         int rx, int ry, int rz) {
  if (lattice == "dimer") {
    // Ignore repetition and use lattice constant as distance.
    return make_unique<dimerBox>(a);
  } else if (lattice == "sc") {
    return make_unique<scBox>(a, rx, ry, rz);
  } else if (lattice == "bcc") {
    return make_unique<bccBox>(a, rx, ry, rz);
  } else if (lattice == "fcc") {
    return make_unique<fccBox>(a, rx, ry, rz);
  } else if (lattice == "diamond") {
    return make_unique<diamondBox>(a, rx, ry, rz);
  } else {
    throw runtime_error("Unknown lattice.");
  }
}


// Helper functions ////////////////////////////////////////////////////

int get_neigh(void *kimmdl, int *mode, int* request,
              int *particle, int *numnei, const int **nei1particle,
              const double **rij) {
  int status;
  unsigned i;
  KIM_API_model& model = **( (KIM_API_model**)kimmdl );
  NeighborListIterator& iter =
    *( (NeighborListIterator*)model.get_data("neighObject", &status) );

  if (*mode == 0) { // Iterator mode
    if (*request == 0) { // Reset
      iter.i = 0;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    } else { // Increment
      i = iter.i;
      if (i >= iter.b.natoms())
        return KIM_STATUS_NEIGH_ITER_PAST_END;
      ++iter.i;
    }
  } else { // Locator mode
    i = *request;
  }

  *particle = i;
  const vector<int>& neighbors = iter.b.get_neigh_list()[i];
  *numnei = neighbors.size();
  *nei1particle = &neighbors[0];
  if (iter.has_rvec)
    *rij = &iter.b.get_neigh_rvec()[i][0];
  else
    *rij = NULL;

  return KIM_STATUS_OK;
}

vector<string> get_free_parameters(KIM_API_model& model) {
  vector<string> retval;
  int n, status;
  char* params = model.get_free_params(&n, &status);
  if (status < KIM_STATUS_OK)
    throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                        + string(" of file ") + string(__FILE__));
  for (int i = 0; i != n; ++i) {
    string s = &params[i*KIM_KEY_STRING_LENGTH];
    retval.push_back(s);
  }
  return retval;
}


// This class contains the box and it also handles KIM interaction. ////
class Compute {
public:
  Compute(const string& lattice, double a,
          bool neigh_iterator, bool neigh_locator, KIMNeigh neighmode)
    : box_(make_box(lattice, a, 1, 1, 1)),//8, 16, 32)),
      iter(*box_, neighmode == KIM_neigh_rvec_f),
      natoms(box_->natoms()),
      ntypes(-1), // Will be set later in the constructor.
      forces(natoms, 3),
      particleEnergy(natoms),
      virial(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      particleVirial(natoms, 6),
      counter(0)
  {
    int status;
    clock_t start_time, stop_time;

    // Assemble KIM descriptor file.
    string descriptor = KIM_DESC;
    if (!neigh_iterator && !neigh_locator)
      throw runtime_error("At least one of neigh_iterator and neigh_locator "
                          "must be true!");
    else if (neigh_iterator && neigh_locator)
      descriptor += "Neigh_BothAccess  flag\n";
    else if (neigh_iterator)
      descriptor += "Neigh_IterAccess  flag\n";
    else if (neigh_locator)
      descriptor += "Neigh_LocaAccess  flag\n";
    switch (neighmode) {
    case KIM_cluster:
      descriptor += "CLUSTER           flag\n";
      break;
    case KIM_mi_opbc_f:
      descriptor += "MI_OPBC_F         flag\n";
      break;
    case KIM_neigh_pure_f:
      descriptor += "NEIGH_PURE_F      flag\n";
      break;
    case KIM_neigh_rvec_f:
      descriptor += "NEIGH_RVEC_F      flag\n";
      break;
    default:
      throw runtime_error("Unsupported neighbor list mode.");
    }

    // Query and store particle type codes.  For this we first get a
    // dummy model instance that can be queried for all the
    // information in the KIM descriptor.  Also add species to the
    // test descriptor.
    KIM_API_model* query;
    status = KIM_API_model_info(&query, MODELNAME); // TODO: hardcoded   
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
    //status = KIM_API_init(&model, TESTNAME, MODELNAME);
    status = KIM_API_string_init(&model, descriptor.c_str(),
                                 MODELNAME); // TODO: hardcoded   
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));

    // Set references to data for KIM.
    model->setm_data(&status, 10*4,
                     "numberOfParticles", 1, &box_->natoms(), 1,
                     "numberParticleTypes", 1, &ntypes, 1,
                     "cutoff", 1, &cutoff, 1,
                     "energy", 1, &energy, 1,
                     "forces", 3*natoms, &forces(0,0), 1,
                     "particleEnergy", natoms, &particleEnergy[0], 1,
                     "virial", 6, &virial(0), 1,
                     "particleVirial", 6*natoms, &particleVirial(0,0), 1,
                     "neighObject", 1, &iter, 1,
                     "boxSideLengths", 3, &box_->box_size()[0], 1
                     );
    if (status < KIM_STATUS_OK) {
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    }

    // Pass methods to KIM.
    status = model->set_method("get_neigh", 1, (func_ptr) &get_neigh);
    if (status < KIM_STATUS_OK) {
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    }

    // Init KIM model.
    status = model->model_init();
    if (status < KIM_STATUS_OK) {
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    }

    // Build neighbor list.
    start_time = clock();
    if (neighmode != KIM_cluster)
      box_->init_neighbor_list(cutoff);
    stop_time = clock();
    neighbor_time = double(stop_time - start_time) / CLOCKS_PER_SEC;

    // Give ghost coordinates.
    const int natoms_all = box_->natoms_all();
    model->setm_data(&status, 2*4, //vvvvvv TODO: this is wrong!   
                     "coordinates", 3*natoms_all, &box_->ghost_coordinates()(0,0), 1,
                     "particleTypes", natoms_all, &box_->ghost_types()[0], 1
                     );
    if (status < KIM_STATUS_OK) {
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));
    }


    // Store parameter names.
    free_parameters = get_free_parameters(*model);
  }

  ~Compute() {
    int status;
    status = model->model_destroy();
    if (status < KIM_STATUS_OK)
      cout << string("KIM error in line ") + to_string(__LINE__)
        + string(" of file ") + string(__FILE__) << endl;
    // TODO: this segfaults...     
    KIM_API_free(&model, &status);
    if (status < KIM_STATUS_OK)
      cout << string("KIM error in line ") + to_string(__LINE__)
        + string(" of file ") + string(__FILE__) << endl;
  }

  void compute() {
    model->model_compute();
  }

  Box& box() {
    return *box_;
  }

  double get_energy() const {
    return energy;
  }

  double get_energy_per_atom() const {
    return energy / natoms;
  }

  double get_neighbor_time() const {
    return neighbor_time;
  }

  /// The stress from virial without kinetic energy contribution in GPa.
  Voigt6<double> get_stress() const {
    double V = box_->get_volume();
    return 160.2177 * virial / V; // GPa
  }

  double get_cutoff() const {
    return cutoff;
  }

  void set_param(const string& name, const string& elem1,
                 const string& elem2, const string& elem3, double value) {
    int status;
    double* p = (double*)model->get_data(name.c_str(), &status); // TODO:     
                                                                 // type      
                                                                 // of        
                                                                 // parameter 
                                                                 // is        
                                                                 // hardcoded.
    if (status < KIM_STATUS_OK)
      throw runtime_error(string("KIM error in line ") + to_string(__LINE__)
                          + string(" of file ") + string(__FILE__));

    map<const string,int>::const_iterator e1 = partcl_type_codes.find(elem1);
    if (e1 == partcl_type_codes.end())
      throw runtime_error(string("unknown species: ") + elem1);

    map<const string,int>::const_iterator e2 = partcl_type_codes.find(elem2);
    if (e2 == partcl_type_codes.end())
      throw runtime_error(string("unknown species: ") + elem2);

    map<const string,int>::const_iterator e3 = partcl_type_codes.find(elem3);
    if (e3 == partcl_type_codes.end())
      throw runtime_error(string("unknown species: ") + elem3);

    Array3D<double> param(p, ntypes, ntypes, ntypes);  // TODO: shape   
                                                       // of parameter  
                                                       // is hardcoded  

    param(e1->second,e2->second,e3->second) = value;

    model->model_reinit();
  }

  //! Objective function for fitting box extents.
  static double objective_box(unsigned n, const double* x, double* grad,
                              void* f_data) {
    Compute& c = *(static_cast<Compute*>(f_data));
    // Scale box.  TODO: We assume the neighbor lists don't change.   
    if (n != 3)
      throw runtime_error("This objective function only works with "
                          "3 parameters.");
    c.box().scale_to(x[0], x[1], x[2]);
    c.compute();
    ++c.counter;
    // Fill gradient.
    if (grad)
      throw runtime_error("Gradient not supported for box optimization.");
    if (false && c.counter % 1 == 0) {
      cout << c.counter << "  " << c.energy << endl;
      c.box().write_to(string("dump.") + to_string(c.counter));
    }
    return c.energy;
  }

  //! Minimize box shape, no constraints.  TODO: non-ortho
  double minimize_box(double ftol_abs, int maxeval) {
    // The box extends are our parameters.
    vector<double> size = {
      box_->box_size()[0], box_->box_size()[1], box_->box_size()[2]
    };
    vector<double> lb = { 0.0, 0.0, 0.0 };
    // TODO: We use a gradient-free algorithm for now, although the
    // gradient should be obtainable from the virial (is it even
    // exactly the virial?).
    double obj_val;
    nlopt::opt optimizer(nlopt::LN_SBPLX, 3);
    optimizer.set_min_objective(Compute::objective_box, this);
    optimizer.set_lower_bounds(lb); // It is a box size, may only be
                                    // positive.
    optimizer.set_initial_step(0.05); // This may not be too big!  
    counter = 0;
    optimizer.set_maxeval(maxeval);
    optimizer.set_ftol_abs(ftol_abs);
    optimizer.optimize(size, obj_val);
    cout << obj_val << " in " << counter << " steps." << endl;
    counter = 0;
    return obj_val;
  }

  //! Objective function for fitting atomic positions.
  static double objective(unsigned n, const double* x, double* grad,
                          void* f_data) {
    Compute& c = *(static_cast<Compute*>(f_data));
    // Write back all positions.  TODO: We assume the neighbor lists
    // don't change.  
    double* coords = &c.box().coordinates()(0,0);
    const Vec3D<double>& box_size = c.box().box_size();
    for (unsigned i = 0; i != n; ++i)
      coords[i] = pmod(x[i], box_size[i % 3]); // TODO: we assume orthorhombic boxes!
    c.box_->update_ghosts(c.box_->nshells(c.cutoff));
    c.compute();
    ++c.counter;
    // Fill gradient.
    if (grad) {
      const double* f = &c.forces(0,0);
      for (unsigned i = 0; i != n; ++i)
        grad[i] = -f[i];
    }
    if (false && c.counter % 1 == 0) {
      cout << c.counter << "  " << c.energy << endl;
      c.box().write_to(string("dump.") + to_string(c.counter));
    }
    return c.energy;
  }

  //! Minimize atomic positions.
  double minimize_positions(double ftol_abs, int maxeval) {
    // Get coordinates as vector.
    Array2D<double>& coords = box_->coordinates();
    const int nparams = coords.extent(0)*coords.extent(1);
    const double* p = &coords(0,0);
    vector<double> coord_vec(p, p+nparams);
    // Other good algorithms:
    //   * NLOPT_LD_LBFGS
    //   * NLOPT_LD_TNEWTON_PRECOND_RESTART
    //   * NLOPT_LD_VAR2
    double obj_val;
    nlopt::opt optimizer(nlopt::LD_TNEWTON_PRECOND_RESTART, nparams);
    optimizer.set_min_objective(Compute::objective, this);
    counter = 0;
    optimizer.set_maxeval(maxeval);
    optimizer.set_ftol_abs(ftol_abs);
    optimizer.optimize(coord_vec, obj_val);
    cout << obj_val << " in " << counter << " steps." << endl;
    counter = 0;
    return obj_val;
  }

  double minimize_positions(double ftol_abs) {
    return minimize_positions(ftol_abs, 10000);
  }

  double minimize_positions() {
    return minimize_positions(0.0001);
  }

private:
  unique_ptr<Box> box_;
  NeighborListIterator iter;
  int natoms, ntypes;
  double cutoff, energy;
  Array2D<double> forces;
  vector<double> particleEnergy;
  Voigt6<double> virial;
  Array2D<double> particleVirial;
  KIM_API_model* model;
  double neighbor_time;
  vector<string> free_parameters; /// Names of the free parameters.
  map<const string,int> partcl_type_codes;
  unsigned counter;
};


// Helpers for main. ///////////////////////////////////////////////////
void print_structure(const string& lattice, double a, KIMNeigh neighmode) {
  printf("%-8s ", lattice.c_str()); fflush(stdout);
  Compute c(lattice, a, true, true, neighmode);
  if (lattice == "fcc") {
    c.set_param("PARAM_FREE_Rc", "Si", "Si", "Si", 2.9);
    c.set_param("PARAM_FREE_Dc", "Si", "Si", "Si", 0.15);
  }
  const clock_t start_time = clock();
  c.compute();
  const clock_t stop_time = clock();
  const Voigt6<double> stress = c.get_stress();
  printf("%8.5f eV/atom       %+.2e  %+.2e  %+.2e   %7d atoms in %6.3f s + %6.3f s\n",
         c.get_energy_per_atom(),
         stress.xx, stress.yy, stress.zz,
         c.box().natoms(),
         c.get_neighbor_time(),
         double(stop_time - start_time) / CLOCKS_PER_SEC);
}



// Main ////////////////////////////////////////////////////////////////

int main() {
  Compute diamond("diamond", 5.429, true, true, KIM_mi_opbc_f);
  //Compute diamond("sc", 2.525, true, true, KIM_mi_opbc_f);
  /*
  diamond.box().coordinates()(0,0) += 0.5;
  diamond.box().coordinates()(10,2) -= 0.7;
  diamond.box().coordinates()(30,1) += 0.2;
  diamond.minimize_positions(0.0001, 10000);
  */
  //diamond.minimize_box(0.00001, 10000);
  cout << "Lattice constants:" << endl;
  cout << " a = " << diamond.box().box_size()[0] / 1 << endl;
  cout << " b = " << diamond.box().box_size()[1] / 1 << endl;
  cout << " c = " << diamond.box().box_size()[2] / 1 << endl;

  diamond.box().init_neighbor_list(diamond.get_cutoff());
  cout << "cutoff = " << diamond.get_cutoff() << endl;
  const vector< vector<int> >& neigh_list = diamond.box().get_neigh_list();
  for (unsigned i = 0; i != diamond.box().natoms(); ++i) {
    cout << "Neighbor distances of atom " << i << ":";
    map<string,int> hist;
    for (int j : neigh_list[i]) {
      const double d = diamond.box().dist(i,j);
      char key[6];
      snprintf(key, 6, "%5.2f", d);
      string k(key);
      hist[k] += 1;
    }
    for (const auto& kv : hist)
      cout << " {" << kv.first << ": " << kv.second << " }";
    cout << endl;
  }

  diamond.compute();
  cout << diamond.get_energy_per_atom() << " eV/atom" << endl;

  /*
  typedef pair<KIMNeigh,string> np;
  vector<np> neigh_modes = {
    //np(KIM_cluster, "cluster"),
    np(KIM_mi_opbc_f, "MI_OPBC"),
    np(KIM_neigh_pure_f, "neighbor list"),
    np(KIM_neigh_rvec_f, "neighbor list with distance vectors")
  };
  for (const pair<KIMNeigh,string>& neigh_mode : neigh_modes) {
    cout << "* Neighbor mode: " << neigh_mode.second << "\n" << endl;
    cout << "struc    E                      pxx        pyy        pzz" << endl;
    print_structure("dimer",   2.232, neigh_mode.first);
    print_structure("sc",      2.525, neigh_mode.first);
    print_structure("bcc",     3.043, neigh_mode.first);
    print_structure("fcc",     3.940, neigh_mode.first);
    print_structure("diamond", 5.429, neigh_mode.first);
    cout << "\n";
  }
  */

  /*
  ofstream outfile;
  outfile.open("diamond.bulkmod.dat");
  outfile << "# V             E           factor\n"
          << "# (ang/atom)    (eV/atom)\n";
  cout << "Calculating rigid bulk modulus for diamond Si..." << flush;
  const clock_t start_time = clock();
  Compute diamond("diamond", 5.429, true, true, KIM_mi_opbc_f);
  const Vec3D<double> size0 = diamond.box().box_size();
  //int j = 0;
  const double deviation = 0.15; // min = 1 - deviation, max = 1 + deviation
  const int halfsteps = 10; // halfsteps*2 + 1 = real steps
  for (int i = -halfsteps; i <= halfsteps; ++i) {
    const double factor = 1.0 + (deviation/halfsteps)*i;
    const Vec3D<double> new_size = factor * size0;
    diamond.box().scale_to(new_size[0], new_size[1], new_size[2]);
    diamond.compute();
    outfile << diamond.box().get_volume() / diamond.box().natoms()
            << " "
            << diamond.get_energy_per_atom()
            << " "
            << factor
            << endl;
    //const string dumpname = string("dump.") + to_string(j++);
    //diamond.box().write_to(dumpname);
  }
  const clock_t stop_time = clock();
  cout << "done in "
       << double(stop_time - start_time) / CLOCKS_PER_SEC
       <<" s, written to diamond.bulkmod.dat" << endl;
  */

  return 0;
}
