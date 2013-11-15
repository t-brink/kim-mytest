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

#ifndef MYTESTCORE_HPP
#define MYTESTCORE_HPP

#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <memory>

#include <KIM_API.h>

#include "ndarray.hpp"

namespace mytest {
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
        @param neighmode The neighbor list mode.
        @param name A name for the box.
    */
    Box(const Vec3D<double>& a, const Vec3D<double>& b, const Vec3D<double>& c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        std::unique_ptr< Array2D<double> > coordinates,
        std::unique_ptr< Array1D<int> > types,
        KIMNeigh neighmode, const std::string& name);

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
        @param neighmode The neighbor list mode.
        @param name A name for the box.

        @todo: lattices with c != a.
    */
    Box(const std::string& lattice, double lattice_const, bool cubic,
        unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        const std::vector<int>& types,
        KIMNeigh neighmode, const std::string& name);

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

        With a MI_OPBC_F neighbor mode, this will also enforce
        periodic boundary conditions but only on the ghost positions
        array, i.e. invisible from the @c positions field.

        Call after manipulating atomic positions even if ghost atoms
        are disabled.
    */
    void update_ghosts();

    /** Update distance vectors.

        Call this after changing the atom positions when using
        NEIGH_RVEC_F. This implicitly calls update_ghosts().

        Calling this method is cheaper than re-calcuting the neighbor
        list with update_neighbor_list().
    */
    void update_ghost_rvecs();

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
     */
    const Vec3D<double>& a;

    /** Constant reference to box vector b.

        This is a three-component vector of a box vector, data is
        continous.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!
     */
    const Vec3D<double>& b;

    /** Constant reference to box vector c.

        This is a three-component vector of a box vector, data is
        continous.

        @todo Actually make sure the internal storage of Vec3D is
        continuous!!  Use array internally, make x/y/z references!
     */
    const Vec3D<double>& c;

    /** Periodic boundary conditions.

        Constant three-component vector that defines if the box is
        periodic in ±a, ±b, and ±c direction.
    */
    const Vec3D<bool>& periodic;

    /** The KIM neighbor mode supported by this box. */
    const KIMNeigh kim_neighbor_mode;

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
    const double* get_neighbor_rvecs_ptr(unsigned i) const {
      return &neigh_rvec_[i][0];
    };

    /** Return const pointer to position (including ghosts).

        This is different from the @c positions field insofar as this
        the data that must be used for actual calculations. It is
        updated by update_neighbor_list() and update_ghosts().  The
        returned pointer may lose validity when calling
        update_neighbor_list().  It is also the position data that
        includes ghost atom positions.

        In short: If you want to change positions use the @c positions
        field, if you want to calculate based on positions use this
        method.
    */
    const double* get_positions_ptr() const {
      return &(*ghost_positions)(0,0);
    }

    /** Return const pointer to types (including ghosts).

        Same arguments as for get_positions() vs @c positions apply.

        @see get_positions()
    */
    const int* get_types_ptr() const {
      return &(*ghost_types)(0);
    }

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

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param factor The factor that the box lengths will be
                      multiplied with.
    */
    void scale(double factor);

    /** Scale box in a, b, and c direction.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param factor_a Multiply box vector a by this factor.
        @param factor_b Multiply box vector b by this factor.
        @param factor_c Multiply box vector c by this factor.
    */
    void scale(double factor_a, double factor_b, double factor_c) {
      deform(Voigt6<double>(factor_a, factor_b, factor_c, 0, 0, 0));
    }

    /** Set box vector lengths.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param len_a Scale box vector a to have this length.
        @param len_b Scale box vector b to have this length.
        @param len_c Scale box vector c to have this length.
    */
    void scale_to(double len_a, double len_b, double len_c) {
      scale(len_a / box_side_lengths_[0],
            len_b / box_side_lengths_[1],
            len_c / box_side_lengths_[2]);
    }


    /** Deform box by given deformation matrix.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        A deformation matrix is defined is the strain tensor plus the
        identity matrix: ε<sub>ij</sub> + δ<sub>ij</sub>.

        @param defmatrix The deformation matrix in Voigt notation.

        @todo Take care that the box stays orthorhombic if
              kim_neighbor_mode is KIM_mi_opbc_f.
    */
    void deform(Voigt6<double> defmatrix);

    /** Deform box to fit the given box vectors.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param new_a The new box vector a.
        @param new_b The new box vector b.
        @param new_c The new box vector c.

        @todo Take care that the box stays orthorhombic if
              kim_neighbor_mode is KIM_mi_opbc_f.
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

  private:
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
    std::vector< std::vector< Vec3D<int> > > neigh_rvec_shell_;
  };


  /** Interface to KIM model free parameters. */
  class FreeParam {
  public:
    FreeParam(const std::string& name,
              const std::vector<unsigned>& shape,
              int kim_index)
      : name(name), rank(shape.size()), size(mult_(shape)), shape(shape),
        kim_index(kim_index)
    {}

    const std::string name;
    const unsigned rank;
    const unsigned size;
    const std::vector<unsigned> shape;
    const int kim_index;

  private:
    /** Helper to get the product of all vector entries. */
    static unsigned mult_(const std::vector<unsigned>& shape) {
      unsigned product = 1;
      for (unsigned factor : shape)
        product *= factor;
      return product;
    }
  };


  /** Main interface class.

      This class contains a simulation box and interacts with KIM.
  */
  class Compute {
  public:
    /** General constructor.

        Take a simulation box and a neighbor mode and prepare KIM.

        @param box The simulation box. This class takes possession of
                   the pointer.
        @param modelname KIM model identifier.
     */
    Compute(std::unique_ptr<Box> box, const std::string& modelname);

    /** Destructor.

        Will deallocate all KIM memory, too.
    */
    ~Compute();

    /** Update the neighbor list and re-register any changed arrays to
        KIM.

        Call after atom positions have changed more than the cutoff
        skin or if the cutoff has changed (increased).
    */
    void update_neighbor_list();

    /** Get box vector a. */
    Vec3D<double> get_a() const { return box_->a; }

    /** Get box vector b. */
    Vec3D<double> get_b() const { return box_->b; }

    /** Get box vector c. */
    Vec3D<double> get_c() const { return box_->c; }

    /** Get code for a particle type.

        @param name String representation (chemical symbol).
        @return Type code.
    */
    int get_particle_type_code(const std::string& name) {
      auto it = partcl_type_codes.find(name);
      if (it != partcl_type_codes.end())
        return it->second;
      else
        throw std::runtime_error("unknown type code");
    }

    /** Get name for a particle type.

        @param code Integer represting that particle type.
        @return Type name.
    */
    std::string get_particle_type_name(int code) {
      auto it = partcl_type_names.find(code);
      if (it != partcl_type_names.end())
        return it->second;
      else
        throw std::runtime_error("unknown type name");
    }

    /** Get position of an atom.

        @param i The index of the atom.
        @return The position of atom i.
    */
    Vec3D<double> get_position(unsigned i) const {
      const double* pos = box_->get_positions_ptr();
      return Vec3D<double>(pos[i*3 + 0], pos[i*3 + 1], pos[i*3 + 2]);
    }

    /** Set position of an atom.

        Does not update neighbor lists.

        @param i The index of the atom.
        @param new_pos The new position of atom i.
        @param update_ghosts Update the actual positions and ghost
                             atom positions of the box. This must be
                             set to true to actually register your
                             change. It may be set to false when
                             updating the positions of several atoms
                             in a row. Then only the last call needs
                             to set this to true. Optional, default is
                             @c true.
    */
    void set_position(unsigned i, const Vec3D<double> new_pos,
                      bool update_ghosts = true){
      box_->positions(i,0) = new_pos[0];
      box_->positions(i,1) = new_pos[1];
      box_->positions(i,2) = new_pos[2];
      if (update_ghosts)
        box_->update_ghost_rvecs();
    }

    /** Set position of an atom.

        Does not update neighbor lists.

        @param i The index of the atom.
        @param x The new x position of atom i.
        @param y The new y position of atom i.
        @param z The new z position of atom i.
        @param update_ghosts Update the actual positions and ghost
                             atom positions of the box. This must be
                             set to true to actually register your
                             change. It may be set to false when
                             updating the positions of several atoms
                             in a row. Then only the last call needs
                             to set this to true. Optional, default is
                             @c true.
    */
    void set_position(unsigned i, double x, double y, double z,
                      bool update_ghosts = true){
      box_->positions(i,0) = x;
      box_->positions(i,1) = y;
      box_->positions(i,2) = z;
      if (update_ghosts)
        box_->update_ghost_rvecs();
    }

    /** Move atom by specified offset.

        Does not update neighbor lists.

        @param i The index of the atom.
        @param offset The offset. Will be added to the position of
                      atom i.
        @param update_ghosts Update the actual positions and ghost
                             atom positions of the box. This must be
                             set to true to actually register your
                             change. It may be set to false when
                             updating the positions of several atoms
                             in a row. Then only the last call needs
                             to set this to true. Optional, default is
                             @c true.
    */
    void move_atom(unsigned i, const Vec3D<double> offset,
                   bool update_ghosts = true){
      box_->positions(i,0) += offset[0];
      box_->positions(i,1) += offset[1];
      box_->positions(i,2) += offset[2];
      if (update_ghosts)
        box_->update_ghost_rvecs();
    }

    /** Move atom by specified offset.

        Does not update neighbor lists.

        @param i The index of the atom.
        @param offset_x The x offset. Will be added to the position of
                        atom i.
        @param offset_y The y offset. Will be added to the position of
                        atom i.
        @param offset_z The z offset. Will be added to the position of
                        atom i.
        @param update_ghosts Update the actual positions and ghost
                             atom positions of the box. This must be
                             set to true to actually register your
                             change. It may be set to false when
                             updating the positions of several atoms
                             in a row. Then only the last call needs
                             to set this to true. Optional, default is
                             @c true.
    */
    void move_atom(unsigned i,
                   double offset_x, double offset_y, double offset_z,
                   bool update_ghosts = true){
      box_->positions(i,0) += offset_x;
      box_->positions(i,1) += offset_y;
      box_->positions(i,2) += offset_z;
      if (update_ghosts)
        box_->update_ghost_rvecs();
    }

    /** Compute values using KIM. */
    void compute();

    /** Get the computed potential energy of the box. */
    double get_energy() const {
      return energy;
    };

    /** Get the computed potential energy per atom of the box. */
    double get_energy_per_atom() const {
      return energy / box_->natoms;
    };

    /** Get the computed virial tensor. */
    Voigt6<double> get_virial() const {
      return virial;
    }

    /** Change model parameter.

        @param param_name KIM name of the free parameter.
        @param indices Indices indicating which element of the
                       parameter array to access (must have zero
                       length for parameters that are scalars).
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.

        @todo There is no way to know the *type* of the free parameter 
        using the KIM API.  We don't check for now.                    
    */
    void set_parameter(const std::string& param_name,
                       const std::vector<unsigned>& indices,
                       double value, bool reinit = true);

    /** Change model parameter.

        @param param_name KIM name of the free parameter.
        @param indices Indices indicating which element of the
                       parameter array to access (must have zero
                       length for parameters that are scalars).
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.

        @todo There is no way to know the *type* of the free parameter 
        using the KIM API.  We don't check for now.                    
    */
    void set_parameter(const std::string& param_name,
                       const std::vector<unsigned>& indices,
                       int value, bool reinit = true);

    /** Optimize box vector lengths to minimize energy.

        Does not update neighbor lists.

        @param ftol_abs Convergence criterium: maximum change of
                        objective function value between two
                        iterations.
        @param maxeval Maximum number of iterations.

        @return The potential energy of the optimized box.
    */
    double optimize_box(double ftol_abs, unsigned maxeval);

    /** Optimize atomic positions to minimize energy.

        Does not update neighbor lists.

        @param ftol_abs Convergence criterium: maximum change of
                        objective function value between two
                        iterations.
        @param maxeval Maximum number of iterations.

        @return The potential energy of the optimized box.
    */
    double optimize_positions(double ftol_abs, unsigned maxeval);

    /** Return number of iterations in last optimization. */
    unsigned get_optimization_steps() const {
      return fit_counter;
    }

    /** Switch the contained Box objects between two Compute objects.

        Both boxes must have the same KIM neighbor list type and use
        the same KIM model.

        @param other Another Compute object with which to switch boxes.
    */
    void switch_boxes(Compute& other);

    /** Store a new box in the object and return the old one.

        Both boxes must have the same KIM neighbor list type.

        This will implicitly call Box::update_neighbor_list() to
        ensure that correct cutoffs are set etc.

        @param new_box The new box.
        @return The old box.
    */
    std::unique_ptr<Box> change_box(std::unique_ptr<Box> new_box);

  private:
    std::unique_ptr<Box> box_;
    const std::string modelname_;
    KIM_API_model* model;
    KIMNeigh neighbor_mode;
    // Model input.
    int ntypes;
    // Model output (constant).
    std::map<std::string,int> partcl_type_codes;
    std::map<int,std::string> partcl_type_names;
    double cutoff;
    std::vector<FreeParam> free_parameters; // TODO: accessors   
    std::map<std::string,FreeParam> free_parameter_map;
    // Computation results, fixed length.
    double energy;
    Voigt6<double> virial;
    // Computation results, length as function of number of atoms +
    // ghost atoms.
    std::vector<double> forces;
    std::vector<double> particleEnergy;
    std::vector<double> particleVirial;
    // Store current index of KIM iteration. This (at least) makes
    // this class not thread safe!
    unsigned kim_iter_pos;
    bool kim_wants_rvec;
    // Some counter used when fitting to get number of iterations (not
    // thread safe obviously).
    unsigned fit_counter;

    /** Get neighbor list function for KIM.

        @param[in] kimmdl The KIM model.
        @param[in] mode Iterator or locator mode.
        @param[in] request The requested index or reset/increment for
                           iterator mode.
        @param[out] particle The actual central atom index.
        @param[out] numnei The number of neighbors.
        @param[out] nei1particle Array of neighbor indices.
        @param[out] rij Array of neighbor distance vectors.

        @return Status.
     */
    static int get_neigh(KIM_API_model** kimmdl,
                         const int* mode, const int* request,
                         int* particle, int* numnei, int** nei1particle,
                         double** rij);

    /** Objective function for optimizing the box.

        This optimizes the length of the box vectors.
    */
    static double obj_func_box(const std::vector<double>& x,
                               std::vector<double>& grad,
                               void* f_data);

    /** Objective function for optimizing the atom positions. */
    static double obj_func_pos(const std::vector<double>& x,
                               std::vector<double>& grad,
                               void* f_data);

    /** set_parameter() implementation. */
    template<typename T>
    void set_parameter_impl(const std::string& param_name,
                            const std::vector<unsigned>& indices,
                            T value, bool reinit);

    /** What it says on the tin. */
    void update_kim_after_box_change();
  };

}

#endif /* MYTESTCORE_HPP */
