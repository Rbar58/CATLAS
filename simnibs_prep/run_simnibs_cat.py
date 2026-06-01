#!/usr/bin/env python3
"""
run_simnibs_cat.py
==================

Run a simple tDCS (transcranial direct current stimulation) electric-field
simulation on the cat head model (cat_head.msh).

>>> RUN THIS ON YOUR OWN COMPUTER, inside the SimNIBS python environment <<<
   (e.g.  `simnibs_python run_simnibs_cat.py`)
The cloud session that generated the mesh can't reach your local SimNIBS.

The mesh `cat_head.msh` is already built and delivered to you, so there is NO
meshmesh step to run. Tissue labels follow SimNIBS' standard numbers
(2=brain[GM] 3=CSF 4=bone 5=skin), so default conductivities apply automatically.

The electrode coordinates below are two well-separated points on the scalp,
auto-computed from this specific model (world/mm coordinates; SimNIBS projects
them onto the nearest scalp surface). Replace them with anatomically meaningful
positions for your study -- they're just a sensible default to get a field.
"""

from simnibs import sim_struct, run_simnibs

# --- electrode positions on the scalp, in mm (world coords of THIS model) -----
ANODE   = [-16.48, 61.19, -46.09]
CATHODE = [ 29.22, 23.29,  32.79]
# extra scalp points you can swap in:
# ALT1 = [-32.24, 31.26, 15.03];  ALT2 = [43.60, 38.33, -31.51]


def main():
    S = sim_struct.SESSION()
    S.subpath = None                 # custom mesh -> no m2m_ folder
    S.fnamehead = "cat_head.msh"
    S.pathfem = "simu_cat_tdcs"      # output folder

    tdcs = S.add_tdcslist()
    tdcs.currents = [1e-3, -1e-3]    # +1 mA anode, -1 mA cathode

    anode = tdcs.add_electrode()
    anode.channelnr = 1
    anode.centre = ANODE
    anode.shape = "ellipse"
    anode.dimensions = [10, 10]      # mm (small -- this is a cat)
    anode.thickness = 2

    cathode = tdcs.add_electrode()
    cathode.channelnr = 2
    cathode.centre = CATHODE
    cathode.shape = "ellipse"
    cathode.dimensions = [10, 10]
    cathode.thickness = 2

    run_simnibs(S)
    print("Done. Open the .msh in simu_cat_tdcs/ with gmsh to see the E-field (normE).")


if __name__ == "__main__":
    main()
