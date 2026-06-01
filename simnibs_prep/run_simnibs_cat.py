#!/usr/bin/env python3
"""
run_simnibs_cat.py
==================

Run a simple tDCS (transcranial direct current stimulation) electric-field
simulation on the synthesized cat head model.

>>> RUN THIS ON YOUR OWN COMPUTER, inside the SimNIBS python environment <<<
   (e.g.  `simnibs_python run_simnibs_cat.py`)
This cannot run in the cloud session that generated the files -- SimNIBS is a
desktop install on your machine.

PREREQUISITE
------------
First turn the labelled volume into a SimNIBS mesh (run once, in a terminal):

    meshmesh head_labels.nii.gz cat_head.msh --voxsize_meshing 0.5

`head_labels.nii.gz` already uses SimNIBS' standard tissue numbers
(1=WM 2=GM 3=CSF 4=bone 5=skin), so default conductivities apply automatically.

The coordinates below are in the SAME world/mm space as head_labels.nii.gz
(computed from the scalp surface). Tweak them to move the electrodes.
"""

from simnibs import sim_struct, run_simnibs

# --- electrode positions on the scalp, in mm (world coords of the model) -----
SCALP_TOP   = [  2.04,  -6.76, 22.12]   # top of head
SCALP_FRONT = [ -1.79,  24.76,  4.99]   # front of head
SCALP_LEFT  = [-22.84,  -7.57,  0.33]
SCALP_RIGHT = [ 22.86,  -8.81, -1.20]


def main():
    S = sim_struct.SESSION()
    S.subpath = None                 # custom mesh -> no m2m_ folder
    S.fnamehead = "cat_head.msh"     # produced by meshmesh
    S.pathfem = "simu_cat_tdcs"      # output folder

    tdcs = S.add_tdcslist()
    tdcs.currents = [1e-3, -1e-3]    # +1 mA anode, -1 mA cathode

    anode = tdcs.add_electrode()
    anode.channelnr = 1
    anode.centre = SCALP_TOP
    anode.shape = "ellipse"
    anode.dimensions = [10, 10]      # mm (small -- this is a cat)
    anode.thickness = 2

    cathode = tdcs.add_electrode()
    cathode.channelnr = 2
    cathode.centre = SCALP_FRONT
    cathode.shape = "ellipse"
    cathode.dimensions = [10, 10]
    cathode.thickness = 2

    run_simnibs(S)
    print("Done. Open the .msh in simu_cat_tdcs/ with gmsh to see the E-field.")


if __name__ == "__main__":
    main()
