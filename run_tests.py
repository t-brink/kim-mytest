#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import itertools
import tempfile
import glob
import shutil
import time

import ase, ase.data, ase.io
import numpy

from mmlib.formats import lammps


# TODO: parse arguments
#
# TODO: extract Tersoff parameters from model and put them in LAMMPS
#       format automatically

elements = ["Si", "C"]
elem_map = {e: i
            for i, e in enumerate(elements, 1)}
elem_map_reverse = {i: e
                    for i, e in enumerate(elements, 1)}

model = "Tersoff_LAMMPS_Erhart_Albe_CSi__MO_903987585848_002"

lammps_potfile = "Tersoff_SiC_ErhartAlbe.params"

lammps_model = """\
pair_style tersoff
pair_coeff * * {} {}
""".format(lammps_potfile, " ".join(elements))

lammps_cmd = ["/home/brink/bin/lmp_serial"]

lattices = {"diamond": {"Si": "5.429", "C": "3.85"},
            "sc":      {"Si": "2.525", "C": "1.80"},
            "bcc":     {"Si": "3.043", "C": "2.15"},
            "fcc":     {"Si": "3.940", "C": "2.80"},
            "random":  {"Si": "1.8",   "C": "1.5"}}

shear = 0.05

commands = {"numer_forces_deriv": "1 1 1",
            "diff_total_energy_vs_particle_energy": "3 3 3",
            "diff_total_virial_vs_particle_virial": "3 3 3",
            "diff_total_virial_vs_virial_from_forces": "3 3 3"}

LAMMPS_TEMPLATE = """
units metal
boundary {px} {py} {pz}
atom_style atomic

box tilt large
read_data {infile}

{masses}
{potential}

timestep 0.001
thermo_style custom step temp pe ke etotal enthalpy press vol lx ly lz xy xz yz pxx pyy pzz pxy pxz pyz cpu cpuremain
thermo 1

compute stress all stress/atom NULL virial

# Convert virial from bar*A^3 to eV
variable xx atom c_stress[1]*6.241508e-7
variable yy atom c_stress[2]*6.241508e-7
variable zz atom c_stress[3]*6.241508e-7
variable xy atom c_stress[4]*6.241508e-7
variable xz atom c_stress[5]*6.241508e-7
variable yz atom c_stress[6]*6.241508e-7

compute pe all pe/atom

dump MyDump all custom 1 {outfile} id type x y z c_pe fx fy fz v_xx v_yy v_zz v_xy v_xz v_yz
dump_modify MyDump format line "%d %d %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g %30.15g"

fix 1 all nve

run 0
"""


with tempfile.TemporaryDirectory() as tmpdir:
    boxinfos = {}

    ########################################################################
    # Run self-consistency tests, while writing all boxes to a temp dir.   #
    ########################################################################
    print("Running self-consistency tests...", end="\n",
          file=sys.stderr, flush=True)
    start = time.time()

    proc = subprocess.Popen(["./mytest"],
                            stdin=subprocess.PIPE,
                            universal_newlines=True)
    def ex(cmd):
        proc.stdin.write(cmd + "\n")

    # Create random box only once per element.
    rand_boxes = {}
    if "random" in lattices:
        for i in itertools.product(sorted(elements),
                                   (True, False), (True, False), (True, False)):
            elem, pbc_x, pbc_y, pbc_z = i
            h_i = hash(i)
            boxname = "randbox-" + ("p" if h_i >= 0 else "m") + "{:x}".format(abs(h_i))
            rand_boxes[i] = boxname
            min_dist = lattices["random"][elem]
            ex("random_box {} {} {} {} {} NEIGH_PURE {}".format(
                boxname, pbc_x, pbc_y, pbc_z, min_dist, elem
            ))

    firstrun = True
    combinations = list(itertools.product(sorted(commands),
                                          sorted(lattices),
                                          sorted(elements),
                                          # PBC
                                          (True, False),
                                          (True, False),
                                          (True, False),
                                          # Cubic        del atom?
                                          (True, False), (True, False),
                                          (-shear, 0, shear),
                                          (-shear, 0, shear),
                                          (-shear, 0, shear)))
    niter = len(combinations)
    last_batch = ()
    for i, combi in enumerate(combinations):
        (command, lattice, elem, pbc_x, pbc_y, pbc_z, cubic,
         del_atom, shear_yz, shear_xz, shear_xy) = combi
        repeat = commands[command]
        latconst = lattices[lattice][elem]
        new_batch = (command, lattice, elem)
        if new_batch != last_batch:
            print("    {:5.1f}%  --  {} {} {}".format(100*i/niter, *new_batch),
                  file=sys.stderr)
            last_batch = new_batch
        # Create the box.
        h_i = hash(combi)
        boxname = "box-" + ("p" if h_i >= 0 else "m") + "{:x}".format(abs(h_i))
        if lattice == "random":
            randboxname = rand_boxes[(elem, pbc_x, pbc_y, pbc_z)]
            ex("copy_box {} {}".format(randboxname, boxname))
        else:
            ex("box {} {} {} {} {} {} {} {} NEIGH_PURE {}".format(
                boxname, lattice, latconst, cubic, repeat,
                pbc_x, pbc_y, pbc_z, elem
            ))
        # Init compute.
        if firstrun:
            ex("model comp {} {}".format(boxname, model))
            firstrun = False
        # Modifications to the box.
        if del_atom:
            ex("delete_atom {}".format(boxname))
        if shear_yz or shear_xz or shear_xy:
            ex("deform_box comp 1 1 1 {} {} {}".format(shear_yz,
                                                       shear_xz,
                                                       shear_xy))
        # Compute.
        p = lambda pbc: "P" if pbc else "O"
        boxinfo = ("{:2s} {:7s} {:1s} {:1s} {:1s} {:7s} "
                   "1 1 1 {:+5.2f} {:+5.2f} {:+5.2f}{}"
                   "".format(elem, lattice,
                             p(pbc_x), p(pbc_y), p(pbc_z),
                             ("cubic" if cubic else "minimal"),
                             shear_yz, shear_xz, shear_xy,
                             (" delete random atom" if del_atom else "")))
        ex("println ==> {:<45s} ".format(command) + boxinfo)
        ex("change_box comp {}".format(boxname))
        # Write box to tmpdir (only for one command, the boxes are the
        # same for all commands; use a replicated one)
        if command == "diff_total_energy_vs_particle_energy":
            boxinfos[boxname] = boxinfo
            pbc_str = p(pbc_x) + p(pbc_y) + p(pbc_z)
            ex(
                "write_box {} {}".format(
                    boxname,
                    os.path.join(tmpdir,
                                 boxname + "-" + pbc_str + ".xyz")
                )
            )
        # Compute.
        if command == "numer_forces_deriv" and lattice == "random":
            # This is too slow, so use the fast (but inexact) version
            command = "numer_forces_deriv_fast"
        ex("{} comp".format(command))

    proc.stdin.close()
    proc.wait()

    del ex
    del proc

    stop = time.time()
    print("done in {:.2f} s.".format(stop-start), file=sys.stderr)

    ########################################################################
    # Run LAMMPS on the boxes.                                             #
    ########################################################################
    print("Running LAMMPS on boxes...", end=" ", file=sys.stderr, flush=True)
    start = time.time()

    empty_boxes = 0
    for box in glob.iglob(os.path.join(tmpdir, "box-*.xyz")):
        fname = os.path.basename(box)
        # Prepare LAMMPS input.
        pbc_str = fname.rsplit(".", 1)[0].split("-")[-1]
        pbc_bool = [{"O": 0, "P": 1}[c]
                    for c in pbc_str]
        pbc = [{"O": "f", "P": "p"}[c]
               for c in pbc_str]
        fname_data = fname.rsplit(".", 1)[0] + ".data"
        outname = fname.rsplit(".", 1)[0] + ".dump"
        masses = "\n".join(
            "mass {} {} # {}".format(
                i,
                ase.data.atomic_masses[ase.data.atomic_numbers[e]],
                e
            )
            for i, e in enumerate(elements, 1)
        )
        lmp_in = LAMMPS_TEMPLATE.format(
            px=pbc[0], py=pbc[1], pz=pbc[2],
            infile=fname_data,
            outfile=outname,
            masses=masses,
            potential=lammps_model
        )
        with open(os.path.join(tmpdir, "lmp.in"), "w") as f:
            f.write(lmp_in)
        # Convert infile.
        ase_box = ase.io.read(box, format="extxyz")
        if not len(ase_box):
            empty_boxes += 1
            continue
        # Ensure that atoms in nonperiodic directions lie inside
        # the box by shifting them a bit and expanding the box in
        # that direction.
        if not all(pbc_bool): # at least one nonperiodic.
            offset = numpy.array([0.0, 0.0, 0.0])
            a = ase_box.cell[0]
            if not pbc_bool[0]:
                offset += 0.25*a
                a *= 1.5
            b = ase_box.cell[1]
            if not pbc_bool[1]:
                offset += 0.25*b
                b *= 1.5
            c = ase_box.cell[2]
            if not pbc_bool[2]:
                offset += 0.25*c
                c *= 1.5
            ase_box.positions += offset
            ase_box.set_cell([a,b,c], scale_atoms=False)
        # Write LAMMPS dump.
        lammps.write_data(os.path.join(tmpdir, fname_data), ase_box,
                          elem_map, ntypes_all=True)
        shutil.copy(lammps_potfile, tmpdir)
        # Run LAMMPS.
        subprocess.run(lammps_cmd + ["-in", "lmp.in"],
                       cwd=tmpdir)
        # Read output and write test file.
        try:
            ase_box, cols = lammps.read_dump(os.path.join(tmpdir, outname),
                                             elem_map_reverse,
                                             extra_columns=True)
        except Exception as e:
            print(file=sys.stderr)
            print("Exception occured when trying to read LAMMPS output:",
                  file=sys.stderr)
            print(e, file=sys.stderr)
            print("See {}".format(tmpdir), file=sys.stderr)
            os._exit(0)
        testfile = fname.rsplit(".", 1)[0] + ".test"
        with open(os.path.join(tmpdir, testfile), "w") as f:
            # Comment line
            f.write(boxinfos[fname.rsplit(".", 1)[0].rsplit("-", 1)[0]] + "\n")

            f.write(str(len(ase_box)) + "\n") # natoms

            f.write(str(ase_box.cell[0,0]) + " ") # a
            f.write(str(ase_box.cell[0,1]) + " ")
            f.write(str(ase_box.cell[0,2]) + " ")
            f.write(str(pbc_bool[0]) + " ")

            f.write(str(ase_box.cell[1,0]) + " ") # b
            f.write(str(ase_box.cell[1,1]) + " ")
            f.write(str(ase_box.cell[1,2]) + " ")
            f.write(str(pbc_bool[1]) + " ")

            f.write(str(ase_box.cell[2,0]) + " ") #c
            f.write(str(ase_box.cell[2,1]) + " ")
            f.write(str(ase_box.cell[2,2]) + " ")
            f.write(str(pbc_bool[2]) + "\n")

            for items in zip(ase_box,
                             cols["c_pe"],
                             cols["fx"],
                             cols["fy"],
                             cols["fz"],
                             cols["v_xx"],
                             cols["v_yy"],
                             cols["v_zz"],
                             cols["v_yz"],  # yz
                             cols["v_xz"],  # xz
                             cols["v_xy"]): # xy
                atom = items[0]
                f.write(atom.symbol + " " +
                        str(atom.position[0]) + " " +
                        str(atom.position[1]) + " " +
                        str(atom.position[2]) + " ")
                f.write(" ".join(str(i) for i in items[1:]))
                f.write("\n")

    stop = time.time()
    print("done in {:.2f} s. Had {} empty boxes."
          "".format(stop-start, empty_boxes),
          file=sys.stderr)

    ########################################################################
    # Test the boxes.                                                      #
    ########################################################################
    print("Comparing models results to LAMMPS results...", end=" ",
          file=sys.stderr, flush=True)
    start = time.time()

    proc = subprocess.Popen(["./mytest"],
                            stdin=subprocess.PIPE,
                            universal_newlines=True)
    def ex(cmd):
        proc.stdin.write(cmd + "\n")

    # Perform init.
    ex("box INITBOX fcc 1.0 true 1 1 1 true true true NEIGH_PURE {}"
       "".format(elements[0]))
    ex("model comp INITBOX {}".format(model))

    # Go.
    for box in glob.iglob(os.path.join(tmpdir, "box-*.test")):
        fname = os.path.basename(box)
        boxname = fname.rsplit(".", 1)[0]
        ex("println")
        ex("check_testfile comp {} {}".format(box, boxname))

    proc.stdin.close()
    proc.wait()

    stop = time.time()
    print("done in {:.2f} s.".format(stop-start), file=sys.stderr)
