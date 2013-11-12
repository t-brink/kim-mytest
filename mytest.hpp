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

#ifndef MYTEST_HPP
#define MYTEST_HPP

#include "ndarray.hpp"

namespace mytest {
  enum KIMNeigh {
    KIM_cluster, KIM_mi_opbc_f, KIM_neigh_pure_f, KIM_neigh_rvec_f
  };

  class Box {
  protected:
    /** Constructor only to be called by subclasses.

      @param n Number of atoms.
      @param a First vector spanning the simulation box.
      @param b Second vector spanning the simulation box.
      @param c Third vector spanning the simulation box.
      @param periodic_a Periodic boundary condition in a-direction.
      @param periodic_b Periodic boundary condition in b-direction.
      @param periodic_c Periodic boundary condition in c-direction.
      @param neigh_list Support neighbor list?
      @param ghosts Prepare ghost atoms? Ignored if @p neigh_list is false.
      @param name A name for the box.
     */
    Box(int n,
        const Vec3D<double>& a, const Vec3D<double>& b, const Vec3D<double>& c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        bool neigh_list, bool ghosts,
        const std::string& name);

  public:
    virtual ~Box() {}

    /** Update ghost atom positions from central atom positions.

        Call after manipulating atomic positions even if ghost atoms
        are disabled.
    */
    void update_ghosts();

    /** Initialize and/or update neighbor list.

        This will allocate memory for ghost atoms if needed.

        @param cutoff The cutoff.
        @param skin Fraction of cutoff that should be added to the
                    cutoff to avoid re-calculating the neighbor list
                    when atoms move small distances.
    */
    void update_neighbor_list(double cutoff, double skin);

    /** Return number of atoms.

        @return Number of atoms without ghost atoms.

        @todo Make this a const reference.
    */
    unsigned natoms() const;

    /** Return number of ghost atoms.

        @return Number of ghost atoms, zero if ghost atoms are disabled.

        @todo Make this a const reference.
    */
    unsigned nghosts() const;

    /** Return number of all atoms, including ghost atoms.

        @return Number of all atoms.

        @todo Make this a const reference.
    */
    unsigned nall() const;

    /** Constant reference to box side lengths.

        This is a three-component vector of side lengths, data is
        continous.

        Be warned that these are just the norms of the box
        vectors. Only for an orthorhombic box do they define the box
        completely.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!

        @todo Do I really want to return a reference? KIM needs one
        but perhaps another solution is better?
    */
    const Vec3D<double>& box_side_lengths;

    /** Constant reference to box vector a.

        This is a three-component vector of a box vector, data is
        continous.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!

        @todo Do I really want to return a reference? KIM needs one
        but perhaps another solution is better?
     */
    const Vec3D<double>& a;

    /** Constant reference to box vector b.

        This is a three-component vector of a box vector, data is
        continous.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!

        @todo Do I really want to return a reference? KIM needs one
        but perhaps another solution is better?
     */
    const Vec3D<double>& b;

    /** Constant reference to box vector c.

        This is a three-component vector of a box vector, data is
        continous.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!

        @todo Do I really want to return a reference? KIM needs one
        but perhaps another solution is better?
     */
    const Vec3D<double>& c;

    /** Atom positions.

        Users are free to manipulate the data in this array.  Changes
        will not be applied until calling update_ghosts().
    */
    Array2D<double> positions;

    /** Calculate volume of the simulation cell.

        @return Volume of the box.
    */
    double calc_volume() const;

    /** Calculate distance between to atoms i and j.

        @param i First atom.
        @param j Second atom, both may be ghost atoms.
        @return Distance.
    */
    double calc_dist(int i, int j) const {
      // Just a convenience wrapper.
      double dx, dy, dz;
      calc_dist(i,j,dx,dy,dz);
    };

    /** Calculate distance between to atoms i and j.

        @param i First atom.
        @param j Second atom, both may be ghost atoms.
        @param[out] dx Will be set to x<sub>j</sub> - x<sub>i</sub>.
        @param[out] dy Will be set to y<sub>j</sub> - y<sub>i</sub>.
        @param[out] dz Will be set to z<sub>j</sub> - z<sub>i</sub>.
        @return Distance.
    */
    double calc_dist(int i, int j, double& dx, double& dy, double& dz) const;

    /** Return const reference to the neighbor list of an atom.

        @param i The index of the central atom, must not be a ghost atom.
        @return A reference to atom i's neighbor list.
    */
    const vector<unsigned>& get_neighbors(unsigned i) const;

    /** Return const pointer to the neighbor list of an atom.

        @param i The index of the central atom, must not be a ghost atom.
        @return A pointer to atom i's neighbor list.
    */
    const unsigned* get_neighbors_ptr(unsigned i) const;

    /** Return const reference to the distance vectors of the
        neighbors of an atom.

        @param i The index of the central atom, must not be a ghost atom.
        @return A reference to a 2D array with distance vectors. In
                the same order as the output of get_neighbors().

        @todo Implementation could be hairy if I want to return a
              reference.  Perhaps the pointer is enough?
    */
    const Array2D<double>& get_neighbor_rvecs(unsigned i) const;

    /** Return const pointer to the distance vectors of the
        neighbors of an atom.

        @param i The index of the central atom, must not be a ghost atom.
        @return A pointer to a 2D array with distance vectors. In
                the same order as the output of get_neighbors().
    */
    const double* get_neighbor_rvecs_ptr(unsigned i) const;

  };

}

#endif /* MYTEST_HPP */
