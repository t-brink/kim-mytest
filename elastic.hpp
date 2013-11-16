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

#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include "mytestcore.hpp"

namespace mytest {

  /** Parameters for the Birch-Murnaghan equation of state. */
  struct BMParams {
    double E0, V0, B0, dB0_dp;
  };

  /** Calculate bulk modulus using an energy-volume curve.

      The bulk modulus is obtained by fitting the Birch-Murnaghan
      equation of state to the energy-volume curve.

      The box that is passed should already be at equilibrium.

      @param compute The Compute object for which to calculate the
                     bulk modulus.
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

      @todo: OpenMP
  */
  BMParams bulk_modulus_energy(Compute& compute,
                               std::vector<double>& volumes,
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

      The box that is passed should already be at equilibrium.

      @param compute The Compute object for which to calculate the
                     bulk modulus.
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

      @todo: OpenMP
  */
  inline
  BMParams bulk_modulus_energy(Compute& compute,
                               double max_strain = 0.05,
                               bool c_to_a = false,
                               bool b_to_a = false,
                               bool positions = false,
                               bool angle_ab = false,
                               bool angle_ac = false,
                               bool angle_bc = false) {
    std::vector<double> volumes;
    std::vector<double> energies;
    return bulk_modulus_energy(compute, volumes, energies, max_strain,
                               c_to_a, b_to_a, positions,
                               angle_ab, angle_ac, angle_bc);
  }

  /** Calculate bulk modulus using a pressure-volume curve.

      The bulk modulus is obtained by fitting the Birch-Murnaghan
      equation of state to the pressure-volume curve.

      The box that is passed should already be at equilibrium.

      @param compute The Compute object for which to calculate the
                     bulk modulus.
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
      @return A Birch-Murnaghan result object. The energy @c E0 will
              be set to zero because it is not contained in the p(V)
              equation of state.

      @todo: OpenMP
  */
  BMParams bulk_modulus_pressure(Compute& compute,
                                 std::vector<double>& volumes,
                                 std::vector<double>& pressures,
                                 double max_strain = 0.05,
                                 bool c_to_a = false,
                                 bool b_to_a = false,
                                 bool positions = false,
                                 bool angle_ab = false,
                                 bool angle_ac = false,
                                 bool angle_bc = false);

  /** Calculate bulk modulus using a pressure-volume curve.

      The bulk modulus is obtained by fitting the Birch-Murnaghan
      equation of state to the pressure-volume curve.

      The box that is passed should already be at equilibrium.

      @param compute The Compute object for which to calculate the
                     bulk modulus.
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
      @return A Birch-Murnaghan result object. The energy @c E0 will
              be set to zero because it is not contained in the p(V)
              equation of state.
  */
  inline
  BMParams bulk_modulus_pressure(Compute& compute,
                                 double max_strain = 0.05,
                                 bool c_to_a = false,
                                 bool b_to_a = false,
                                 bool positions = false,
                                 bool angle_ab = false,
                                 bool angle_ac = false,
                                 bool angle_bc = false) {
    std::vector<double> volumes;
    std::vector<double> pressures;
    return bulk_modulus_pressure(compute, volumes, pressures, max_strain,
                                 c_to_a, b_to_a, positions,
                                 angle_ab, angle_ac, angle_bc);
  }

}


#endif /* ELASTIC_HPP */
