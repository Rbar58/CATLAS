#!/usr/bin/env python3
"""
build_head_segmentation.py
==========================

Turn the brain-only CATLAS data into an *approximate full-head* tissue
segmentation that SimNIBS can mesh with `meshmesh`.

WHY THIS EXISTS
---------------
SimNIBS simulates how TMS / tES electric fields spread through a head. The
field is dominated by the skull (resistive) and scalp (conductive), which sit
*outside* the brain. CATLAS only contains brain tissue (grey matter, white
matter, CSF) -- it has no skull and no scalp. So before SimNIBS can do anything
useful we have to invent approximate skull and scalp layers around the brain.

WHAT IT DOES
------------
1. Reads CATLAS's tissue probability maps (TPM.nii):
      channel 0 = grey matter, 1 = white matter, 2 = CSF, 3 = non-brain
2. Builds a hard brain segmentation (each voxel -> WM/GM/CSF by largest prob).
3. Grows concentric shells outward from the brain:
      + a thin extra CSF margin (subarachnoid space)
      + a skull shell
      + a scalp shell
4. Writes everything into ONE labelled NIfTI using SimNIBS' STANDARD tissue
   numbers, so SimNIBS applies its default conductivities automatically:

      1 = White matter
      2 = Grey matter
      3 = CSF
      4 = Bone (skull)
      5 = Scalp (skin)

   Output: head_labels.nii.gz   (0.5 mm, ready for `meshmesh`)

IMPORTANT CAVEAT
----------------
The skull and scalp here are IDEALISED SHELLS, not a real cat skull. Absolute
field magnitudes will be approximate. This model is for learning the SimNIBS
pipeline and qualitative work, not for publication-grade dosimetry.
"""

import numpy as np
import nibabel as nib
from scipy import ndimage

# ----------------------------------------------------------------------------
# Tunable shell thicknesses (in voxels; CATLAS voxels are 0.5 mm)
# ----------------------------------------------------------------------------
CSF_MARGIN_VOX = 1   # ~0.5 mm extra CSF around the brain
SKULL_VOX      = 3   # ~1.5 mm skull
SCALP_VOX      = 3   # ~1.5 mm scalp
BRAIN_PROB_THRESHOLD = 0.5   # voxel is "brain" if P(GM+WM+CSF) >= this

# SimNIBS standard tissue labels
WM, GM, CSF, BONE, SKIN = 1, 2, 3, 4, 5

def main():
    tpm_img = nib.load("../TPM.nii")
    tpm = tpm_img.get_fdata()          # shape (150,192,96,4)
    gm, wm, csf, nb = tpm[..., 0], tpm[..., 1], tpm[..., 2], tpm[..., 3]

    seg = np.zeros(gm.shape, dtype=np.int16)

    # --- 1. brain interior: pick the most probable tissue per voxel ---------
    brain_p = gm + wm + csf
    brain_mask = brain_p >= BRAIN_PROB_THRESHOLD
    brain_mask = ndimage.binary_fill_holes(brain_mask)

    stack = np.stack([wm, gm, csf], axis=-1)      # index 0->WM,1->GM,2->CSF
    winner = np.argmax(stack, axis=-1)            # 0,1,2
    label_of = np.array([WM, GM, CSF])
    seg[brain_mask] = label_of[winner[brain_mask]]

    # --- 2. grow shells outward --------------------------------------------
    def shell(current_mask, n_vox):
        """Dilate current_mask by n_vox and return the NEW voxels added."""
        grown = ndimage.binary_dilation(current_mask, iterations=n_vox)
        grown = ndimage.binary_fill_holes(grown)
        return grown & ~current_mask, grown

    current = brain_mask

    new, current = shell(current, CSF_MARGIN_VOX)
    seg[new & (seg == 0)] = CSF       # subarachnoid CSF margin

    new, current = shell(current, SKULL_VOX)
    seg[new & (seg == 0)] = BONE      # skull shell

    new, current = shell(current, SCALP_VOX)
    seg[new & (seg == 0)] = SKIN      # scalp shell

    # --- 3. save ------------------------------------------------------------
    out = nib.Nifti1Image(seg, tpm_img.affine, tpm_img.header)
    out.header.set_data_dtype(np.int16)
    nib.save(out, "head_labels.nii.gz")

    # --- 4. report ----------------------------------------------------------
    names = {WM: "white matter", GM: "grey matter", CSF: "CSF",
             BONE: "skull/bone", SKIN: "scalp/skin"}
    vox_ml = np.prod(tpm_img.header.get_zooms()[:3]) / 1000.0  # mm^3 -> mL
    print("Wrote head_labels.nii.gz  shape", seg.shape,
          " voxel", tuple(round(float(z), 2) for z in tpm_img.header.get_zooms()[:3]), "mm")
    print(f"{'label':>5}  {'tissue':<14} {'voxels':>10} {'volume (mL)':>12}")
    for lab in (WM, GM, CSF, BONE, SKIN):
        n = int((seg == lab).sum())
        print(f"{lab:>5}  {names[lab]:<14} {n:>10} {n*vox_ml:>12.2f}")


if __name__ == "__main__":
    main()
