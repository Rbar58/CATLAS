# Cat head model for SimNIBS — starter kit

This folder turns the **brain-only CATLAS** atlas into an **approximate full
head** that [SimNIBS](https://simnibs.github.io) can mesh and simulate
(TMS / tES electric fields).

## ⚠️ Read this first (the honest caveats)

- CATLAS contains **only brain** (grey matter, white matter, CSF). It has **no
  skull and no scalp**. Electric-field models depend heavily on the skull and
  scalp, so those layers had to be **invented** here by growing shells outward
  from the brain.
- The skull/scalp are therefore **idealised shells, not a real cat skull**.
  Use this for **learning the pipeline and qualitative exploration**, not for
  publication-grade dosimetry.
- SimNIBS' automatic segmentation (`charm`) is built for **human** heads and
  will not work on a cat — that's exactly why we build the segmentation
  manually here instead.

## What's in this folder

| File | Purpose |
|------|---------|
| `build_head_segmentation.py` | Builds the labelled head volume from `../TPM.nii`. |
| `run_simnibs_cat.py` | Example tDCS simulation script (run on YOUR machine). |
| `README.md` | This guide. |

Running `build_head_segmentation.py` produces **`head_labels.nii.gz`** — a
0.5 mm labelled head (1=WM, 2=GM, 3=CSF, 4=skull, 5=scalp). It isn't checked
into git (it's a derived binary); you regenerate it in step 3 below. The labels
use **SimNIBS' standard tissue numbers**, so SimNIBS applies its default
conductivities automatically — no extra setup needed.

## Step-by-step (on your own computer, where SimNIBS is installed)

1. **Install SimNIBS** if you haven't: https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html

2. **Get the whole repo on your computer** (this `simnibs_prep` folder *and*
   `TPM.nii` one level up). The cloud session that generated this can't reach
   your local SimNIBS.

3. **Generate the labelled head volume** (open a terminal in this folder):
   ```
   pip install nibabel numpy scipy      # one-time, if you don't have them
   python build_head_segmentation.py
   ```
   This reads `../TPM.nii` and writes `head_labels.nii.gz`.

4. **Build the mesh:**
   ```
   meshmesh head_labels.nii.gz cat_head.msh --voxsize_meshing 0.5
   ```
   This produces `cat_head.msh`, a tetrahedral head mesh. (`--voxsize_meshing
   0.5` matches the 0.5 mm data and resolves the thin cat-skull shell.)

5. **Look at the mesh** to sanity-check it:
   ```
   gmsh cat_head.msh
   ```
   You should see nested brain → skull → scalp surfaces.

6. **Run the example simulation:**
   ```
   simnibs_python run_simnibs_cat.py
   ```
   Results land in `simu_cat_tdcs/`. Open the resulting `.msh` in `gmsh` to see
   the electric field (`normE`) on the brain.

## Adjusting things

- **Move the electrodes / coil:** edit the coordinates near the top of
  `run_simnibs_cat.py`. They're in the same mm space as `head_labels.nii.gz`.
  Pre-computed scalp landmarks (top / front / left / right) are included.
- **TMS instead of tDCS:** swap `add_tdcslist()` for `add_tmslist()` and set a
  coil `.fnamecoil` plus a position matrix — see the SimNIBS docs.
- **Thicker / thinner skull or scalp:** change `SKULL_VOX` / `SCALP_VOX` at the
  top of `build_head_segmentation.py` and re-run it (needs `nibabel`, `numpy`,
  `scipy`).
- **Custom mesh limitation:** because there's no `m2m_` folder, EEG-style
  positions ("Cz") and MNI mapping aren't available — give electrode/coil
  placements as **coordinates**, as the script does.

## References

- SimNIBS custom meshing: https://simnibs.github.io/simnibs/build/html/tutorial/advanced/custom_meshing.html
- SimNIBS head meshing overview: https://simnibs.github.io/simnibs/build/html/tutorial/head_meshing.html
- CATLAS source: Stolzberg et al. (2017), *J. Comp. Neurology* — cat brain atlas.
