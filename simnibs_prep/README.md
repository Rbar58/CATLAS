# Cat head model for SimNIBS

A subject-specific cat head model for **SimNIBS** electric-field simulation
(TMS / tES), built from a **real whole-head MRI** anchored to the **CATLAS**
brain atlas. The tetrahedral mesh was generated in-pipeline (see
`realmri_pipeline.py`); you only need to run the *simulation* on your machine.

> There are two pipelines in this folder:
> - **`realmri_pipeline.py`** — *recommended.* Uses a real whole-head MRI +
>   CATLAS. This produced the delivered `cat_head.msh`.
> - `build_head_segmentation.py` — older fallback that synthesises a head from
>   the brain-only atlas *without* any real MRI. Keep only if you have no scan.

## ⚠️ Honest quality notes (first-pass model)

This is a **proof-of-concept / first-pass** model — good for validating the
pipeline and qualitative work, not yet publication-grade dosimetry:

- **Skull is approximate.** T1 MRI cannot image bone, so the skull is a
  geometric shell grown around the brain at a fixed ~1.5 mm thickness. For
  accurate dosimetry you need a **CT** (or a dual-echo/UTE scan) to place real
  bone. The skull dominates tES/TMS fields, so treat absolute values with care.
- **Single brain compartment.** This in-vivo T1 has negligible grey/white
  contrast (cluster means ~250 vs ~262), so GM and WM are merged into one
  `brain` label rather than fabricating a split. Add WM later from a better
  scan or by warping the CATLAS WM prior.
- **Brain mask** is registration-derived (CATLAS→subject). The subject scan's
  header had its A-P and S-I axes effectively swapped, which rotated the warped
  brain ~90°; this is now corrected by an **anatomical header fix** before
  registration (see "Orientation" below) and the boundary refined with a
  cross-correlation (`SyNCC`) deformable (brain ≈ 28.7 mL). Still QC every new
  scan — the axis mapping was determined for *this* subject.
- **Other tissues follow the brain.** Scalp comes from the real MRI; CSF is a
  GMM on the subject's N4 intensities *inside* the brain; the bone shell is
  positioned *around* the brain. So they are only as well-placed as the brain —
  which is why fixing the brain orientation was the prerequisite for the rest.
- SimNIBS' automatic segmentation (`charm`) is **human-only** and will not work
  on a cat — that is why the segmentation is built explicitly here.

### Orientation

The subject affine is rewritten so the voxel axes carry their true anatomical
meaning — axis0 = L-R, axis1 = A-P (anterior +Y), axis2 = S-I (superior +Z) —
*before* registration, so SyN has no gross rotation to fight. This was verified
against the anatomy in three orthogonal views (sagittal: anterior up; coronal:
ventral/dorsal across; horizontal/axial: eyes visible, anterior to one side).
See `anatomical_reorient()` in `realmri_pipeline.py`.

## Delivered files (sent to you directly; not committed — they derive from your unpublished MRI)

| File | What it is |
|------|------------|
| `cat_head.msh` | The tetrahedral head mesh (SimNIBS). Open in `gmsh`. |
| `cat_tissues.nii.gz` | The labelled volume it was meshed from (2=brain 3=CSF 4=bone 5=scalp). |
| `seg_qc.png`, `seg_qc_contour.png` | QC: tissue maps / boundaries on the MRI. |

Labels use **SimNIBS standard tissue numbers**, so default conductivities apply
automatically (note: WM=1 is absent; brain is 2=GM with GM conductivity).

## Run a simulation (on your computer, where SimNIBS is installed)

1. Put `cat_head.msh` and `run_simnibs_cat.py` in the same folder.
2. Edit the electrode coordinates at the top of `run_simnibs_cat.py` if you want
   (they default to two well-separated scalp points auto-computed from the model;
   they're world/mm coordinates and get projected onto the scalp).
3. Run:
   ```
   simnibs_python run_simnibs_cat.py
   ```
   Results land in `simu_cat_tdcs/`. Open that `.msh` in `gmsh` to see the
   electric field (`normE`) on the brain.

Because this is a custom mesh (no `m2m_` folder), give electrode/coil positions
as **coordinates**, not EEG names like "Cz".

## Reproduce or refine the model

You normally don't need to — the mesh is delivered — but the full pipeline is
here and reproducible:

```bash
# one-time: enable SimNIBS meshing in a headless Linux/py3.11 env
bash install_simnibs_meshing.sh          # see script header for what it does

# then, with your whole-head MRI:
export LD_LIBRARY_PATH="$(python3 -c 'import sys,os;print(os.path.join(sys.prefix,"lib"))'):$LD_LIBRARY_PATH"
python realmri_pipeline.py --mri your_head.nii --atlas-dir .. --out ./out
# -> out/cat_tissues.nii.gz  and  out/cat_head.msh
```

To improve quality: add a CT for the skull, hand-correct the brain mask in
ITK-SNAP (`out/cat_tissues.nii.gz`), then re-run `--stage 3` to re-mesh.

## References

- SimNIBS custom meshing: https://simnibs.github.io/simnibs/build/html/tutorial/advanced/custom_meshing.html
- SimNIBS head meshing overview: https://simnibs.github.io/simnibs/build/html/tutorial/head_meshing.html
- CATLAS source: Stolzberg et al. (2017), *J. Comp. Neurology* — cat brain atlas.
