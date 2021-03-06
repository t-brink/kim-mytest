/*
  Copyright (c) 2013,2017,2018,2019,2020 Tobias Brink

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
#include <random>

#include "KIM_SimulatorHeaders.hpp"

#include "ndarray.hpp"
#include "utils.hpp"

/** The main namespace. */
namespace mytest {
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
        @param name A name for the box.
    */
    Box(const Vec3D<double>& a, const Vec3D<double>& b, const Vec3D<double>& c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        std::unique_ptr< Array2D<double> > coordinates,
        std::unique_ptr< Array1DInit<std::string> > types,
        const std::string& name);

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
        @param name A name for the box.

        @todo: lattices with c != a.
    */
    Box(const std::string& lattice, double lattice_const, bool cubic,
        unsigned repeat_a, unsigned repeat_b, unsigned repeat_c,
        bool periodic_a, bool periodic_b, bool periodic_c,
        const std::vector<std::string>& types,
        const std::string& name);

    /** Copy constructor that changes the name. */
    Box(const Box& other, const std::string& new_name);

    /** Copy constructor. */
    Box(const Box& other);

    virtual ~Box() {}

    /** Create a randomly filled box

        @todo: document params

     */
    static
    std::unique_ptr<Box> random_box(double a, double b, double c,
                                    bool periodic_a,
                                    bool periodic_b,
                                    bool periodic_c,
                                    double min_dist,
                                    const std::string& atomtype,
                                    const std::string& name,
                                    std::mt19937& rng);

    /** Create a randomly filled box with random types from list

        @todo: document params

     */
    static
    std::unique_ptr<Box> random_box(double a, double b, double c,
                                    bool periodic_a,
                                    bool periodic_b,
                                    bool periodic_c,
                                    double min_dist,
                                    const std::vector<std::string>& atomtypes,
                                    const std::string& name,
                                    std::mt19937& rng);

    /** Initialize and/or update neighbor list.

        This will allocate memory for ghost atoms if needed and
        includes a call to update_ghosts().

        @todo should wrap atoms outside the box back into it!

        @param cutoff The cutoff.
        @param skin Fraction of cutoff that should be added to the
                    cutoff to avoid re-calculating the neighbor list
                    when atoms move small distances.
        @param typemap Mapping of atom type strings to the model's
                       type codes.
        @return If the return value is true, then any previous
                pointers to ghost atom memory are invalid (that means
                re-setting KIM data for positions and types at least).
    */
    bool update_neighbor_list(double cutoff, double skin,
                              const std::map<std::string,int>& typemap);

    /** Update ghost atom positions from central atom positions.

        With a MI_OPBC_F neighbor mode, this will also enforce
        periodic boundary conditions but only on the ghost positions
        array, i.e. invisible from the @c positions field.

        Call after manipulating atomic positions even if ghost atoms
        are disabled.

        @param typemap Mapping of atom type strings to the model's
                       type codes.
    */
    void update_ghosts(const std::map<std::string,int>& typemap);

    /** Return unique_ptr to new box with atom i deleted.

        @param i         The atom that should be deleted.
        @param new_name  Name for the new box.

     */
    std::unique_ptr<Box> delete_atom(unsigned i, const std::string& name) const;

    /** Return unique_ptr to new box with atom i deleted.

        The new box will have the same name as the old one.

        @param i         The atom that should be deleted.

     */
    std::unique_ptr<Box> delete_atom(unsigned i) const {
      return delete_atom(i, name_);
    }

    /** Return unique_ptr to new box, repeated along the given directions.

        @param repeat_a  Repeat this many times in a direction.
        @param repeat_b  Repeat this many times in b direction.
        @param repeat_c  Repeat this many times in c direction.
        @param new_name  Name for the new box.

    */
    std::unique_ptr<Box> repeat(unsigned repeat_a,
                                unsigned repeat_b,
                                unsigned repeat_c,
                                const std::string& name) const;

    /** Return unique_ptr to new box, repeated along the given directions.

        The new box will have the same name as the old one.

        @param repeat_a  Repeat this many times in a direction.
        @param repeat_b  Repeat this many times in b direction.
        @param repeat_c  Repeat this many times in c direction.

    */
    std::unique_ptr<Box> repeat(unsigned repeat_a,
                                unsigned repeat_b,
                                unsigned repeat_c) const {
      return repeat(repeat_a, repeat_b, repeat_c, name_);
    }

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

        Contains atom types as strings.

        Users are free to manipulate the data in this array.  Changes
        will not be applied until calling update_ghosts().  Changing
        the size produces undefined results!
    */
    Array1DInit<std::string> types;

    // Calculate/return values /////////////////////////////////////////

    /** Calculate volume of the simulation cell.

        @return Volume of the box.
    */
    double calc_volume() const {
      return dot(a_, cross(b_, c_));
    }

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

    /** Return const pointer to type codes (including ghosts).

        Same arguments as for get_positions() vs @c positions
        apply. Further difference compared to the @c types field is
        that these are the type codes (int, not string) as the model
        wants them.

        @see get_positions()
    */
    const int* get_types_ptr() const {
      return &(*ghost_types)(0);
    }

    /** Return const pointer to "contributing" array.

        Same arguments as for get_positions() vs @c positions
        apply. Further difference compared to the @c types field is
        that these are the type codes (int, not string) as the model
        wants them.

        @see get_positions()
    */
    const int* get_contributing_ptr() const {
      return &contributing_[0];
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

    /** Write to file in extended XYZ format.

        @param output An output stream.
    */
    void write_to(std::ostream& output) const;

    // Manipulate box //////////////////////////////////////////////////

    /** Uniformly scale box.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param factor The factor that the box lengths will be
                      multiplied with.
        @param typemap Mapping of atom type strings to the model's
                       type codes.
    */
    void scale(double factor,
               const std::map<std::string,int>& typemap) {
      deform(Voigt6<double>(factor, factor, factor, 0, 0, 0), typemap);
    }

    /** Scale box in a, b, and c direction.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param factor_a Multiply box vector a by this factor.
        @param factor_b Multiply box vector b by this factor.
        @param factor_c Multiply box vector c by this factor.
        @param typemap Mapping of atom type strings to the model's
                       type codes.
    */
    void scale(double factor_a, double factor_b, double factor_c,
               const std::map<std::string,int>& typemap);

    /** Set box vector lengths.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        @param len_a Scale box vector a to have this length.
        @param len_b Scale box vector b to have this length.
        @param len_c Scale box vector c to have this length.
        @param typemap Mapping of atom type strings to the model's
                       type codes.
    */
    void scale_to(double len_a, double len_b, double len_c,
                  const std::map<std::string,int>& typemap) {
      scale(len_a / box_side_lengths_[0],
            len_b / box_side_lengths_[1],
            len_c / box_side_lengths_[2],
            typemap);
    }


    /** Deform box by given deformation matrix.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

        A deformation matrix is defined is the strain tensor plus the
        identity matrix: ε<sub>ij</sub> + δ<sub>ij</sub>.

        @param defmatrix The deformation matrix in Voigt notation.
        @param typemap Mapping of atom type strings to the model's
                       type codes.
    */
    void deform(Voigt6<double> defmatrix,
                const std::map<std::string,int>& typemap);

    /** Deform box to fit the given box vectors.

        Includes scaling of atom positions.  Implicitly calls
        update_ghost_rvecs().

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
    std::vector<int> contributing_; // false if ghost
  };


  /** Interface to KIM model free parameters. */
  class FreeParam {
  public:
    FreeParam(const std::string& name,
              const std::string& description,
              KIM::DataType data_type,
              int size,
              int kim_index)
      : name(name), description(description),
        data_type(data_type), size(size), kim_index(kim_index)
    {}

    const std::string name;
    const std::string description;
    const KIM::DataType data_type;
    const unsigned size;
    const int kim_index;
  };


  /** Parameters for the Birch-Murnaghan equation of state. */
  struct BMParams {
    double E0, V0, B0, dB0_dp;
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
    Compute(std::unique_ptr<Box> box, const std::string& modelname,
            KIM::LengthUnit length_unit = KIM::LENGTH_UNIT::A,
            KIM::EnergyUnit energy_unit = KIM::ENERGY_UNIT::eV,
            KIM::ChargeUnit charge_unit = KIM::CHARGE_UNIT::e,
            KIM::TemperatureUnit temperature_unit = KIM::TEMPERATURE_UNIT::K,
            KIM::TimeUnit time_unit = KIM::TIME_UNIT::ps);

    /** Destructor.

        Will deallocate all KIM memory, too.
    */
    ~Compute();

    /** Update the neighbor list and re-register any changed arrays to
        KIM.

        Call after atom positions have changed more than the cutoff
        skin or if the cutoff has changed (increased).

        @param force_ptr_update Force re-registering pointers to KIM.
    */
    void update_neighbor_list(bool force_ptr_update = false);

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

    /** Get number of particle types in the box. */
    unsigned get_number_of_particle_types() {
      return partcl_type_codes.size();
    }

    /** Return number of atoms in the simulation box. */
    unsigned get_natoms() const {
      return box_->natoms;
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
        box_->update_ghosts(partcl_type_codes);
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
        box_->update_ghosts(partcl_type_codes);
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
        box_->update_ghosts(partcl_type_codes);
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
        box_->update_ghosts(partcl_type_codes);
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
    void deform(Voigt6<double> defmatrix) {
      box_->deform(defmatrix, partcl_type_codes);
    }

    /** Compute values using KIM.

        Output values that are not supported by the model are simply
        zero.

        Due to the way periodic boundaries are implemented in CLUSTER
        mode (namely using ghost atoms), we need particleEnergy to
        actually calculate energy and particleVirial to actually
        calculate virial. If a model doesn't provide those energy and
        virial are simply set to zero in that case.
    */
    void compute();

    /** Get the computed potential energy of the box. */
    double get_energy() const {
      return energy;
    };

    /** Get the computed potential energy of atom i. */
    double get_energy(unsigned i) const {
      return particleEnergy[i];
    };

    /** Get the computed potential energy per atom of the box. */
    double get_energy_per_atom() const {
      return energy / box_->natoms;
    };

    /** Get the force of atom i, component dim */
    double get_force(unsigned i, unsigned dim) const {
      if (dim >= 3)
        throw std::runtime_error("invalid dim");
      if (i >= box_->natoms)
        throw std::runtime_error("invalid atom id");
      return forces[3*i + dim];
    }

    /** Get the computed virial tensor. */
    Voigt6<double> get_virial() const {
      return virial;
    }

    /** Get the computed virial tensor for atom i. */
    Voigt6<double> get_virial(unsigned i) const {
      if (i >= box_->natoms)
        throw std::runtime_error("invalid atom id");
      return Voigt6<double>(particleVirial[6*i + 0],
                            particleVirial[6*i + 1],
                            particleVirial[6*i + 2],
                            particleVirial[6*i + 3],
                            particleVirial[6*i + 4],
                            particleVirial[6*i + 5]);
    }

    /** Get the global virial tensor computed from forces.

        This can be used to verify the global virial tensor from the
        model.

     */
    Voigt6<double> get_virial_from_forces() const {
      return Voigt6<double>(global_virial_from_forces_xx,
                            global_virial_from_forces_yy,
                            global_virial_from_forces_zz,
                            global_virial_from_forces_yz,
                            global_virial_from_forces_xz,
                            global_virial_from_forces_xy);
    }

    /** Get the global virial tensor computed using process_dEdr.

        This can be used to verify the global virial tensor from the
        model and process_dEdr.

     */
    Voigt6<double> get_virial_from_dEdr() const {
      return virial_from_dEdr;
    }

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        @param param_name KIM name of the free parameter.
        @param index Index indicating which element of the
                     parameter array to access.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const unsigned index,
                       double value, bool reinit = true);

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        @param param_name KIM name of the free parameter.
        @param index Index indicating which element of the
                     parameter array to access.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const unsigned index,
                       int value, bool reinit = true);

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        This method assumes that the parameter is a N*N row-major
        array, where N is the number of supported species. It also
        assumes that the species code can be used to index into that
        array.

        @param param_name KIM name of the free parameter.
        @param species1 The element name for i.
        @param species2 The element name for j.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const std::string& species1,
                       const std::string& species2,
                       double value, bool reinit = true);

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        This method assumes that the parameter is a N*N*N row-major
        array, where N is the number of supported species. It also
        assumes that the species code can be used to index into that
        array.

        @param param_name KIM name of the free parameter.
        @param species1 The element name for i.
        @param species2 The element name for j.
        @param species3 The element name for k.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const std::string& species1,
                       const std::string& species2,
                       int value, bool reinit = true);

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        This method assumes that the parameter is a N*N*N row-major
        array, where N is the number of supported species. It also
        assumes that the species code can be used to index into that
        array.

        @param param_name KIM name of the free parameter.
        @param species1 The element name for i.
        @param species2 The element name for j.
        @param species3 The element name for k.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const std::string& species1,
                       const std::string& species2,
                       const std::string& species3,
                       double value, bool reinit = true);

    /** Change model parameter.

        Neighbor lists will not be updated automatically.

        This method assumes that the parameter is a N*N*N row-major
        array, where N is the number of supported species. It also
        assumes that the species code can be used to index into that
        array.

        @param param_name KIM name of the free parameter.
        @param species1 The element name for i.
        @param species2 The element name for j.
        @param species3 The element name for k.
        @param value The new value (must match the type of the
                     parameter).
        @param reinit The new parameter may not actually be used or
                      used correctly until this method was called with
                      reinit = true. Setting it to false can be used
                      when updating several parameters in a row to
                      save computation time. Optional, default is @c true.
                      If true, the cutoff will also be updated from KIM.
    */
    void set_parameter(const std::string& param_name,
                       const std::string& species1,
                       const std::string& species2,
                       const std::string& species3,
                       int value, bool reinit = true);


    int get_parameter_int(const std::string& param_name,
                          const unsigned index);

    double get_parameter_double(const std::string& param_name,
                                const unsigned index);

    /** Write model parameters to files.

        @param directory The directory into which the files are
                         written. Has to exist and be writable.
        @param name The name of the new parameterized KIM model.
     */
    void write_kim_parameters(const std::string& directory,
                              const std::string& name) const;

    /** Optimize box vector lengths to minimize energy.

        Does not update neighbor lists.

        @param ftol_abs Convergence criterium: maximum change of
                        objective function value between two
                        iterations. Optional, default is 0.001.
        @param maxeval Maximum number of iterations. Optional, default
                       is 10000.
        @param isotropic If this is true, the box will only be scaled
                         by a single factor, i.e. c/a and b/a are
                         constant. Optional, default is @c false.
        @return The potential energy of the optimized box.
    */
    double optimize_box(double ftol_abs = 0.001,
                        unsigned maxeval = 10000,
                        bool isotropic = false);

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

    /** Calculate bulk modulus using an energy-volume curve.

        The bulk modulus is obtained by fitting the Birch-Murnaghan
        equation of state to the energy-volume curve.

        The simulation box should already be at equilibrium (positions
        and box vectors).

        Previous compute() results will be invalidated.

        @param[out] volumes Will contain volume data points after the call.
        @param[out] energies Will contain energy data points after the call.
        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.05 (5%).
        @param c_to_a Optimize c/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      e.g. cubic structures and costs computing
                      time). Optional, default @c false.
        @param b_to_a Optimize b/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      many structures and costs computing
                      time). Optional, default @c false.
        @param positions Optimize atomic positions. Optional, default
                         @c false.
        @param angle_ab Optimize angle between a and b box
                        vector. Optional, default @c false.
        @param angle_ac Optimize angle between a and c box
                        vector. Optional, default @c false.
        @param angle_bc Optimize angle between b and c box
                        vector. Optional, default @c false.
        @return A Birch-Murnaghan result object.
    */
    BMParams bulk_modulus_energy(std::vector<double>& volumes,
                                 std::vector<double>& energies,
                                 double max_strain = 0.05,
                                 bool c_to_a = false,
                                 bool b_to_a = false,
                                 bool positions = false,
                                 bool angle_ab = false,
                                 bool angle_ac = false,
                                 bool angle_bc = false);

    /** Calculate bulk modulus using an energy-volume curve.

        The bulk modulus is obtained by fitting the Birch-Murnaghan
        equation of state to the energy-volume curve.

        The simulation box should already be at equilibrium (positions
        and box vectors).

        Previous compute() results will be invalidated.

        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.05 (5%).
        @param c_to_a Optimize c/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      e.g. cubic structures and costs computing
                      time). Optional, default @c false.
        @param b_to_a Optimize b/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      many structures and costs computing
                      time). Optional, default @c false.
        @param positions Optimize atomic positions. Optional, default
                         @c false.
        @param angle_ab Optimize angle between a and b box
                        vector. Optional, default @c false.
        @param angle_ac Optimize angle between a and c box
                        vector. Optional, default @c false.
        @param angle_bc Optimize angle between b and c box
                        vector. Optional, default @c false.
        @return A Birch-Murnaghan result object.
    */
    BMParams bulk_modulus_energy(double max_strain = 0.05,
                                 bool c_to_a = false,
                                 bool b_to_a = false,
                                 bool positions = false,
                                 bool angle_ab = false,
                                 bool angle_ac = false,
                                 bool angle_bc = false) {
      std::vector<double> volumes;
      std::vector<double> energies;
      return bulk_modulus_energy(volumes, energies, max_strain,
                                 c_to_a, b_to_a, positions,
                                 angle_ab, angle_ac, angle_bc);
    }

    /** Calculate bulk modulus using an pressure-volume curve.

        The bulk modulus is obtained by fitting the Birch-Murnaghan
        equation of state to the pressure-volume curve.

        The simulation box should already be at equilibrium (positions
        and box vectors).

        Previous compute() results will be invalidated.

        @param[out] volumes Will contain volume data points after the call.
        @param[out] pressures Will contain pressure data points after the call.
        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.05 (5%).
        @param c_to_a Optimize c/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      e.g. cubic structures and costs computing
                      time). Optional, default @c false.
        @param b_to_a Optimize b/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      many structures and costs computing
                      time). Optional, default @c false.
        @param positions Optimize atomic positions. Optional, default
                         @c false.
        @param angle_ab Optimize angle between a and b box
                        vector. Optional, default @c false.
        @param angle_ac Optimize angle between a and c box
                        vector. Optional, default @c false.
        @param angle_bc Optimize angle between b and c box
                        vector. Optional, default @c false.
        @return A Birch-Murnaghan result object.
    */
    BMParams bulk_modulus_pressure(std::vector<double>& volumes,
                                   std::vector<double>& pressures,
                                   double max_strain = 0.05,
                                   bool c_to_a = false,
                                   bool b_to_a = false,
                                   bool positions = false,
                                   bool angle_ab = false,
                                   bool angle_ac = false,
                                   bool angle_bc = false);

    /** Calculate bulk modulus using an pressure-volume curve.

        The bulk modulus is obtained by fitting the Birch-Murnaghan
        equation of state to the pressure-volume curve.

        The simulation box should already be at equilibrium (positions
        and box vectors).

        Previous compute() results will be invalidated.

        @param[out] volumes Will contain volume data points after the call.
        @param[out] pressures Will contain pressure data points after the call.
        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.05 (5%).
        @param c_to_a Optimize c/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      e.g. cubic structures and costs computing
                      time). Optional, default @c false.
        @param b_to_a Optimize b/a ratio for the different volumes (this
                      can be disabled because it is not important in
                      many structures and costs computing
                      time). Optional, default @c false.
        @param positions Optimize atomic positions. Optional, default
                         @c false.
        @param angle_ab Optimize angle between a and b box
                        vector. Optional, default @c false.
        @param angle_ac Optimize angle between a and c box
                        vector. Optional, default @c false.
        @param angle_bc Optimize angle between b and c box
                        vector. Optional, default @c false.
        @return A Birch-Murnaghan result object.
    */
    BMParams bulk_modulus_pressure(double max_strain = 0.05,
                                   bool c_to_a = false,
                                   bool b_to_a = false,
                                   bool positions = false,
                                   bool angle_ab = false,
                                   bool angle_ac = false,
                                   bool angle_bc = false) {
      std::vector<double> volumes;
      std::vector<double> pressures;
      return bulk_modulus_pressure(volumes, pressures, max_strain,
                                   c_to_a, b_to_a, positions,
                                   angle_ab, angle_ac, angle_bc);
    }

    /** Calculate elastic constant c_ij.

        The constant is obtained by fitting to σ_i = c_ij·ε_j.  To
        make this work, the simulation box should already be at
        equilibrium (positions and box vectors).

        This method uses a Voigt notation where
        1 -> xx
        2 -> yy
        3 -> zz
        4 -> yz
        5 -> xz
        6 -> xy,
        which is the same as KIM.

        Previous compute() results will be invalidated.

        @param i First index, must be in [1..6].
        @param j Second index, must be in [1..6].
        @param[out] strain Will contain strain data points after the call.
        @param[out] stress Will contain stress data points after the call.
        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.025 (2.5%).
        @param positions Optimize atomic positions. For some
                         deformation modes this is
                         essential. Optional, default @c false.
    */
    double elastic_constant(unsigned i, unsigned j,
                            std::vector<double>& strain,
                            std::vector<double>& stress,
                            double max_strain = 0.025,
                            bool positions = false);

    /** Calculate elastic constant c_ij.

        The constant is obtained by fitting to σ_i = c_ij·ε_j.  To
        make this work, the simulation box should already be at
        equilibrium (positions and box vectors).

        This method uses a Voigt notation where
        1 -> xx
        2 -> yy
        3 -> zz
        4 -> yz
        5 -> xz
        6 -> xy,
        which is the same as KIM.

        Previous compute() results will be invalidated.

        @param i First index, must be in [1..6].
        @param j Second index, must be in [1..6].
        @param max_strain The maximum strain put on the
                          structure. Optional, default is 0.025 (2.5%).
        @param positions Optimize atomic positions. For some
                         deformation modes this is
                         essential. Optional, default @c false.
    */
    double elastic_constant(unsigned i, unsigned j,
                            double max_strain = 0.025,
                            bool positions = false) {
      std::vector<double> strain;
      std::vector<double> stress;
      return elastic_constant(i, j, strain, stress, max_strain, positions);
    }

    /** Derive force numerically from the energy and return maximum
        deviation from model forces.

        The maximum deviation is for a single dimension of a single
        atom (e.g. atom 12's force in y direction).

        This uses Ridders' method and is slow.
    */
    double max_diff_force();

    /** Derive force numerically from the energy and return maximum
        deviation from model forces.

        The maximum deviation is for a single dimension of a single
        atom (e.g. atom 12's force in y direction).

        This uses the central difference method and is fast (but has a
        larger error).
    */
    double max_diff_force_fast();

    /** Switch the contained Box objects between two Compute objects.

        Will implicitly call update_neighbor_list() for both Compute
        objects.

        @param other Another Compute object with which to switch boxes.
    */
    void switch_boxes(Compute& other);

    /** Store a new box in the object and return the old one.

        Will implicitly call update_neighbor_list().

        @param new_box The new box.
        @return The old box.
    */
    std::unique_ptr<Box> change_box(std::unique_ptr<Box> new_box);

    /** Return a copy of the internal simulation box. */
    std::unique_ptr<Box> copy_box() const {
      return make_unique<Box>(*box_);
    }

    /** Scale the box.

        @see Box::scale
    */
    void scale_box(double factor_a, double factor_b, double factor_c) {
      box_->scale(factor_a, factor_b, factor_c, partcl_type_codes);
    }

    /** Scale the box.

        @see Box::scale
    */
    void scale_box(double factor) {
      box_->scale(factor, partcl_type_codes);
    }

  private:
    std::unique_ptr<Box> box_;
    const std::string modelname_;
    KIM::Model* model;
    KIM::ComputeArguments* compute_arguments;
    int nall_kim; // KIM wants an int for the number of atoms, so we
                  // maintain a copy :-/
    // Inputs and outputs supported by the model.
    bool has_reinit;
    bool has_write_params; // TODO: not handled, yet   
    bool has_energy;
    bool has_particleEnergy;
    bool has_forces;
    bool has_virial;
    bool has_particleVirial;
    bool has_process_dEdr;
    // Model output (constant).
    std::map<std::string,int> partcl_type_codes;
    std::map<int,std::string> partcl_type_names;
    double cutoff;
    std::map<std::string,FreeParam> free_parameter_map;
    // Computation results, fixed length.
    double energy;
    Voigt6<double> virial;
    Voigt6<double> virial_from_dEdr;
    // Computation results, length as function of number of atoms +
    // ghost atoms.
    std::vector<double> forces;
    std::vector<double> particleEnergy;
    std::vector<double> particleVirial;
    // For testing: global virial computed from forces.
    double global_virial_from_forces_xx,
           global_virial_from_forces_yy,
           global_virial_from_forces_zz,
           global_virial_from_forces_yz,
           global_virial_from_forces_xz,
           global_virial_from_forces_xy;
    // Store current index of KIM iteration. This (at least) makes
    // this class not thread safe!
    unsigned kim_iter_pos;
    // Some counter used when fitting to get number of iterations (not
    // thread safe obviously).
    unsigned fit_counter;

    /** Get neighbor list function for KIM.

        @param[in] compute_ptr The Compute object.
        @param[in] number_of_neighlists Ignored for now.  
        @param[in] cutoffs Ignored for now.  
        @param[in] neighbor_list_idx Ignored for now.  
        @param[in] i The central atom index.
        @param[out] n_neighs The number of neighbors.
        @param[out] neighlist Array of neighbor indices.

        @return Status.
     */
    static int get_neigh(void * const compute_ptr,
                         const int number_of_neighlists,
                         const double * const cutoffs,
                         const int neighbor_list_idx,
                         const int i,
                         int * const n_neighs,
                         const int ** const neighlist);

    /** Handle dE/dr. Called from the kim model.

        @param compute_ptr The Compute object
        @param dEdr        Value of dEdr.
        @param r           Pair distance.
        @param vec_r       Pair distance vector.
        @param i           Central atom index.
        @param j           Bonded atom index.

        @return Status.
     */
    static int process_dEdr(void * const compute_ptr,
                            const double dEdr,
                            const double r,
                            const double * const vec_r,
                            const int i,
                            const int j);

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

    /** Find parameter object for given parameter name. */
    const FreeParam& get_param_obj(const std::string& param_name) const {
      const auto it = free_parameter_map.find(param_name);
      if (it == free_parameter_map.end())
        throw std::runtime_error("Unknown free parameter: " + param_name);
      return it->second;
    }

    /** Assuming a N*N array, get flat index for i,j

        @param i Outer index.
        @param j Inner index.
        @param size Size of the flat array.
     */
    static unsigned conv_index(unsigned i, unsigned j, unsigned size) {
      const unsigned n = std::sqrt(size);
      if (n*n != size)
        throw std::runtime_error("\"size\" does not belong to a N*N array.");
      return i*n + j;
    }

    /** Assuming a N*N*N array, get flat index for i,j,k

        @param i Outer index.
        @param j "Middle" index.
        @param k Inner index.
        @param size Size of the flat array.
     */
    static unsigned conv_index(unsigned i, unsigned j, unsigned k,
                               unsigned size) {
      const unsigned n = std::cbrt(size);
      if (n*n*n != size)
        throw std::runtime_error("\"size\" does not belong to a N*N*N array.");
      return i*n*n + j*n + k;
    }

    /** set_parameter() implementation. */
    template<typename T>
    void set_parameter_impl(const FreeParam& param,
                            const unsigned index,
                            T value, bool reinit);

    /** get_parameter() implementation. */
    template<typename T>
    T get_parameter_impl(const std::string& param_name,
                         const unsigned index);

    /** What it says on the tin. */
    void update_kim_after_box_change();
  };

}

#endif /* MYTESTCORE_HPP */
