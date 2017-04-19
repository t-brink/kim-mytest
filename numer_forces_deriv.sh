#!/bin/bash

declare -A lattices

element="Si"
model="Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002"
#lattices=( ["diamond"]="5.429"
#           ["sc"]="2.525"
#           ["bcc"]="3.043"
#           ["fcc"]="3.940" )
lattices=( ["fcc"]="2.22" )
firstrun="yes"

(
  for lattice in "${!lattices[@]}"; do
    latconst="${lattices[$lattice]}"
    for pbcx in "true" "false"; do
      for pbcy in "true" "false"; do
        for pbcz in "true" "false"; do
          boxname="$lattice-$pbcx-$pbcy-$pbcz"
          echo "print === $lattice $latconst $pbcx $pbcy $pbcz ==="
          echo "box $boxname $lattice $latconst true 1 1 1 $pbcx $pbcy $pbcz NEIGH_PURE $element"
          if [[ ! -z $firstrun ]]; then
            echo "model comp $boxname $model"
            firstrun=""
          else
            echo "change_box comp $boxname"
          fi
          echo "numer_forces_deriv comp"
          echo "print"
        done
      done
    done
  done
) | ./mytest
