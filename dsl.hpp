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

#ifndef DSL_HPP
#define DSL_HPP

#include <vector>
#include <map>
#include <string>
#include <memory>

#include "core.hpp"

namespace mytest {
  /** Tokenize, splitting at spaces. */
  std::vector<std::string> tokenize(const std::string& input);

  /** Parse and execute a line of a simple command language.

      The language simply consist of a statement per line. The first
      word is always the subroutine to be executed followed by space
      separated arguments.

      Commands:

      box <name> <lattice> <lattice constant> <cubic?>
          <repeat_a> <repeat_b> <repeat_c>
          <periodic_a?> <periodic_b?> <periodic_c?>
          <neighbor_mode>
          <type1> [<type2> [<type3> [...]]]

         Initialize a Box object.

      model <name> <box> <modelname>

         Initialize a Compute object.

      compute <computer_name>

         Calculates and prints energy and stress (from virial) if
         supported.

      optimize_box <computer_name> [<isotropic?>]

         Optimize box vector lengths and print new lattice vectors, as
         well as energy and stress. Default for isotropic is false.

      optimize_positions <computer_name>

         Optimize atomic positions and print energy and stress.

      bulk_modulus_energy <computer_name> <max_strain>
                          [<c_to_a?> [<b_to_a?> [<positions>
                          [<angle_ab?> [<angle_ac?> [<angle_bc>]]]]]]

         Calculate bulk modulus using energy-volume curves.  Print E0,
         V0, B0, B0'.

      bulk_modulus_pressure <computer_name> <max_strain>
                            [<c_to_a?> [<b_to_a?> [<positions?>
                            [<angle_ab?> [<angle_ac?> [<angle_bc?>]]]]]]

         Calculate bulk modulus using pressure-volume curves.  Print
         V0, B0, B0'.

      stiffness_tensor <computer_name> <max_strain> <positions?>

         Calculate and print all elastic constants c_ij.

      switch_boxes <computer1_name> <computer2_name>

         Switch boxes of both Compute objects.

      change_box <computer_name> <box_name>

         Change box in computer to a new box.

      scale <model_name> <factor_a> [<factor_b> <factor_c>]

         Scale box by given factor(s).

      Arguments with a question mark are booleans, denoted using
      "true"/"false".
  */
  void parse(std::string command,
             std::map< std::string,std::unique_ptr<Box> >& boxes,
             std::map<std::string,Compute>& computes);
}


#endif /* DSL_HPP */
