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

#include <fstream>
#include <vector>
#include <string>
#include <memory>

#include "ndarray.hpp"

namespace mytest {
  /** Utility function that should be in C++14 but not yet in C++11. */
  template<typename T, typename ...Args>
  std::unique_ptr<T> make_unique(Args&& ...args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  /** KIM neighbor mode */
  enum KIMNeigh {
    KIM_cluster, KIM_mi_opbc_f, KIM_neigh_pure_f, KIM_neigh_rvec_f
  };

  /** Class to hold simulation box and neighbor lists.

      @todo Disabled neighbor list leads to segfaults because the
      class does not take into account that there are no neighbor
      lists and tries to access them.

      @todo Way to disable rvec in the constructor.

      @todo Copy constructor, move constructor, and semantics.

      @todo All methods that manipulate the box (scale etc.) should
      also have constant versions that return a new box that has the
      given manipulation done.  This can be useful in some
      circumstances.
  */
  class Box {
  public:
    // Initialization related //////////////////////////////////////////

    /** General constructor.

        This constructor takes unique_ptr arguments to assume
        ownership of position/type arrays. This is done to save memory
        and to avoid copying this data. If you need to copy the data
        you must do that manually beforehand.

        @param a First vector spanning the simulation box.
        @param b Second vector spanning the simulation box.
        @param c Third vector spanning the simulation box.
        @param periodic_a Periodic boundary condition in a-direction.
        @param periodic_b Periodic boundary condition in b-direction.
        @param periodic_c Periodic boundary condition in c-direction.
        @param coordinates Unique pointer to n×3 coordinate array.
        @param types Unique pointer to n component array containing
                     the atom types. Must have the same length as the
                     first dimension of @p coordinates.
        @param neigh_list Support neighbor list?
        @param name A name for the box.
    */
    Box(const Vec3D<double>& a, const Vec3D<double>& b, const Vec3D<double>& c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        std::unique_ptr< Array2D<double> > coordinates,
        std::unique_ptr< Array1D<int> > types,
        bool neigh_list, const std::string& name);

    /** Lattice constructor.

        Construct a given lattice and repeat in specified way.

        @param lattice A string specifying which lattice should be
                       produced. Supported lattices with a basis of 1
                       are @c sc, @c bcc, @c fcc, and @c diamond.
                       Supported lattices with a basis of 2 are @c
                       dimer, @c B1 (NaCl), @c B2 (CsCl), and @c B3
                       (zincblende).
        @param lattice_const The lattice constant a. In case of a
                             dimer this is interpreted as the bond
                             length.
        @param cubic Produce a cubic unit cell instead of the
                     primitive unit cell.
        @param repeat_a Repeat unit cell this many times in a direction.
        @param repeat_b Repeat unit cell this many times in b direction.
        @param repeat_c Repeat unit cell this many times in c direction.
        @param periodic_a Periodic boundary condition in a-direction.
        @param periodic_b Periodic boundary condition in b-direction.
        @param periodic_c Periodic boundary condition in c-direction.
        @param coordinates Unique pointer to n×3 coordinate array.
        @param types The types of the atoms on the lattice. The length
                     of this array must match the basis of the given
                     lattice.
        @param neigh_list Support neighbor list?
        @param name A name for the box.

        @todo: lattices with c != a.
    */
    Box(const std::string& lattice, double lattice_const, bool cubic,
        unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        const std::vector<int>& types,
        bool neigh_list, const std::string& name);

    virtual ~Box() {}

    /** Initialize and/or update neighbor list.

        This will allocate memory for ghost atoms if needed and
        includes a call to update_ghosts().

        @param cutoff The cutoff.
        @param skin Fraction of cutoff that should be added to the
                    cutoff to avoid re-calculating the neighbor list
                    when atoms move small distances.
        @return If the return value is true, then any previous
                pointers to ghost atom memory are invalid (that means
                re-setting KIM data for positions and types at least).
    */
    bool update_neighbor_list(double cutoff, double skin);

    /** Update ghost atom positions from central atom positions.

        Call after manipulating atomic positions even if ghost atoms
        are disabled.
    */
    void update_ghosts();

    // Public data /////////////////////////////////////////////////////

    /** Constant reference to box side lengths.

        This is a three-component vector of side lengths, data is
        continous.

        Be warned that these are just the norms of the box
        vectors. Only for an orthorhombic box do they define the box
        completely.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!
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

    /** Periodic boundary conditions.

        Constant three-component vector that defines if the box is
        periodic in ±a, ±b, and ±c direction.
    */
    const Vec3D<bool>& periodic;

    /** Number of atoms without ghost atoms. */
    const unsigned& natoms;

    /** Number of ghost atoms. */
    const unsigned& nghosts;

    /** Number of all atoms, including ghost atoms. */
    const unsigned& nall;

    /** Atom positions.

        Users are free to manipulate the data in this array.  Changes
        will not be applied until calling update_ghosts().
    */
    Array2D<double> positions;

    /** Atom types.

        Users are free to manipulate the data in this array.  Changes
        will not be applied until calling update_ghosts().  Changing
        the size produces undefined results!
    */
    Array1D<int> types;

    // Calculate/return values /////////////////////////////////////////

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
      return calc_dist(i,j,dx,dy,dz);
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
    const std::vector<int>& get_neighbors(unsigned i) const {
      return neigh_list_[i];
    }

    /** Return const pointer to the neighbor list of an atom.

        @param i The index of the central atom, must not be a ghost atom.
        @return A pointer to atom i's neighbor list.
    */
    const int* get_neighbors_ptr(unsigned i) const {
      return &get_neighbors(i)[0];
    }

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

    // IO //////////////////////////////////////////////////////////////

    /** Write to file.

        @param filename The name of the output file (will be silently
                        overwritten).
    */
    void write_to(const std::string& filename) const {
      std::ofstream f(filename);
      write_to(f);
    };

    /** Write to file.

        @param output An output stream.

        @todo: use some useful format (POSCAR needs sorting by atom type :( )
    */
    void write_to(std::ostream& output) const;

    // Manipulate box //////////////////////////////////////////////////

    /** Uniformly scale box.

        Includes scaling of atom positions.

        @param factor The factor that the box lengths will be
                      multiplied with.
    */
    void scale(double factor);

    /** Scale box in a, b, and c direction.

        Includes scaling of atom positions.

        @param factor_a Multiply box vector a by this factor.
        @param factor_b Multiply box vector b by this factor.
        @param factor_c Multiply box vector c by this factor.
    */
    void scale(double factor_a, double factor_b, double vector_c);

    /** Set box vector lengths.

        Includes scaling of atom positions.

        @param len_a Scale box vector a to have this length.
        @param len_b Scale box vector b to have this length.
        @param len_c Scale box vector c to have this length.
    */
    void scale_to(double len_a, double len_b, double len_c);

    /** Deform box by given deformation matrix.

        Includes scaling of atom positions.

        A deformation matrix is defined is the strain tensor plus the
        identity matrix: ε<sub>ij</sub> + δ<sub>ij</sub>.

        @param defmatrix The deformation matrix in Voigt notation.
    */
    void deform(const Voigt6<double>& defmatrix);

    /** Deform box to fit the given box vectors.

        Includes scaling of atom positions.

        @param new_a The new box vector a.
        @param new_b The new box vector b.
        @param new_c The new box vector c.
    */
    void deform_to(const Vec3D<double>& new_a,
                   const Vec3D<double>& new_b,
                   const Vec3D<double>& new_c);

  protected:
    /** Return number of atoms in the unit cell of a given lattice.

        See the constructor for details of the parameters.
     */
    static unsigned atoms_per_unit_cell(const std::string& lattice,
                                        bool cubic);

    /** Return the number of different species in the unit cell of a
        given lattice.

        See the constructor for details of the parameters.
     */
    static unsigned species_per_unit_cell(const std::string& lattice);

    /** Calculate the number of ghost atom shells needed. */
    Vec3D<unsigned> calc_number_of_ghost_shells(double cutoff) const;

  //private:    
  public:
    Vec3D<double> box_side_lengths_;
    Vec3D<double> a_, b_, c_;
    Vec3D<bool> periodic_;
    unsigned natoms_, nghosts_, nall_;
    std::string name_;
    Vec3D<unsigned> ghost_shells;
    std::unique_ptr< Array2D<double> > ghost_positions;
    std::unique_ptr< Array1D<int> > ghost_types;
    std::vector< std::vector<int> > neigh_list_;
    std::vector< std::vector<double> > neigh_rvec_;
  };

}

#endif /* MYTEST_HPP */
