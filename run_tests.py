#!/usr/bin/env python3

import subprocess
import argparse
import itertools


# TODO: parse arguments

model = "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002"

elements = ["Si", "C"]

lattices = {"diamond": "5.429",
            "sc": "2.525",
            "bcc": "3.043",
            "fcc": "3.940"}

shear = 0.05

commands = {"numer_forces_deriv": "1 1 1",
            "diff_total_energy_vs_particle_energy": "3 3 3",
            "diff_total_virial_vs_particle_virial": "3 3 3",
            "diff_total_virial_vs_virial_from_forces": "3 3 3"}


proc = subprocess.Popen(["./mytest"],
                        stdin=subprocess.PIPE,
                        universal_newlines=True)
def ex(cmd):
    proc.stdin.write(cmd + "\n")

firstrun = True
for i in itertools.product(sorted(commands), sorted(lattices), sorted(elements),
                           (True, False), (True, False), (True, False),
                           (True, False), (True, False),
                           (-shear, 0, shear),
                           (-shear, 0, shear),
                           (-shear, 0, shear)):
    (command, lattice, elem, pbc_x, pbc_y, pbc_z, cubic,
     del_atom, shear_yz, shear_xz, shear_xy) = i
    repeat = commands[command]
    latconst = lattices[lattice]
    # Create the box.
    h_i = hash(i)
    boxname = "box-" + ("p" if h_i >= 0 else "m") + "{:x}".format(abs(h_i))
    ex("box {} {} {} {} {} {} {} {} NEIGH_PURE {}\n".format(
        boxname, lattice, latconst, cubic, repeat, pbc_x, pbc_y, pbc_z, elem
    ))
    # Modifications to the box.
    if del_atom:
        ex("delete_atom {}\n".format(boxname))
    if shear_yz or shear_xz or shear_xy:
        ex("deform_box comp 1 1 1 {} {} {}".format(shear_yz,
                                                   shear_xz,
                                                   shear_xy))
    # Compute.
    p = lambda pbc: "P" if pbc else "O"
    ex(
        "println ==> {:<45s} {:2s} {:7s} {:1s} {:1s} {:1s} {:7s} "
        "1 1 1 {:+5.2f} {:+5.2f} {:+5.2f}{}"
        "\n".format(command, elem, lattice,
                    p(pbc_x), p(pbc_y), p(pbc_z),
                    ("cubic" if cubic else "minimal"),
                    shear_yz, shear_xz, shear_xy,
                    (" delete random atom" if del_atom else ""))
    )
    if firstrun:
        ex("model comp {} {}\n".format(boxname, model))
        firstrun = False
    else:
        ex("change_box comp {}\n".format(boxname))
    ex("{} comp\n".format(command))

proc.stdin.close()
proc.wait()
