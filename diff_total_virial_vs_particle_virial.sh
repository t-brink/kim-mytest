#!/bin/bash

declare -A lattices

element="Si"
model="Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002"
lattices=( ["diamond"]="5.429"
           ["sc"]="2.525"
           ["bcc"]="3.043"
           ["fcc"]="3.940" )
firstrun="yes"
command="diff_total_virial_vs_particle_virial"

(
  for lattice in "${!lattices[@]}"; do
    latconst="${lattices[$lattice]}"
    for pbcx in "true" "false"; do
      for pbcy in "true" "false"; do
        for pbcz in "true" "false"; do
          for cubic in "true" "false"; do
            if [[ $cubic = "true" ]]; then
              c="cubic"
            else
              c="minimal"
            fi
            boxname="$lattice-$pbcx-$pbcy-$pbcz-$c"
            echo "$(
              printf "println %7s %5s - %5s %5s %5s - %7s ::" \
                     $lattice $latconst $pbcx $pbcy $pbcz $c
            )"
            echo "box $boxname $lattice $latconst $cubic 3 3 3 $pbcx $pbcy $pbcz NEIGH_PURE $element"
            if [[ ! -z $firstrun ]]; then
              echo "model comp $boxname $model"
              firstrun=""
            else
              echo "change_box comp $boxname"
            fi
            echo "$command comp"
            # Remove a random atom and re-run.
            echo "delete_atom $boxname"
            echo "$( printf "println %43s :: " "random atom removed" )"
            echo "$command comp"
          done
        done
      done
    done
  done
) | ./mytest