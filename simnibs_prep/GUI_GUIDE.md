# Running the cat model in the SimNIBS GUI (no coordinates needed)

You don't have to touch Python or coordinates — the graphical workflow works
with `cat_head.msh`. Confirmed against SimNIBS 4.6. (Electrode positions are
placed by **clicking on the 3D head**, which is what you're used to.)

## 1. Launch the GUI
In a terminal:
```
simnibs_gui
```
(or use the SimNIBS GUI shortcut the installer created).

## 2. Load the head mesh
In the top of the window there are three fields:

| Field | What to do |
|-------|------------|
| **m2m Folder** | **Leave EMPTY.** (That folder only exists for human `charm` runs; this custom cat mesh doesn't have one and doesn't need one.) |
| **Head Mesh** | Click **Browse** → select `cat_head.msh`. |
| **Output folder** | Click **Browse** → pick where results should go. |

The 3D head (the scalp surface) appears on the right. If it doesn't render,
your machine's OpenGL is the issue, not the mesh — tell me and I'll give you the
scripted route instead.

## 3. Add a stimulation
- For tDCS/tACS: click **"Add tDCS Poslist"** → a tDCS tab appears.
- For TMS: click **"Add TMS Poslist"** instead (coil placement is analogous).

## 4. Place electrodes by clicking on the head (tDCS)
In the tDCS tab:
1. Click **Add** to create an electrode row.
2. **Double-click the electrode's position cell** → an interactive 3D head
   window opens. **Click on the scalp** where you want the electrode centre.
   (For a directional/rectangular electrode it also asks for a y-direction
   point — click a second spot to set orientation.)
3. Set the electrode **shape** (ellipse/rectangle), **size** in mm
   (start small — e.g. 10×10 mm for a cat), **thickness**, and the **current**
   (e.g. +0.001 A on one electrode, −0.001 A on the other; they must sum to 0).
4. Repeat for the second electrode.

> Sizing reality check: this is a cat head (~5 cm), so human-sized 5×5 cm pads
> won't fit. Keep electrodes small.

## 5. Run
Click **Run**. SimNIBS builds the electrodes into the mesh and solves the field.

## 6. Look at the result
When it finishes it opens (or you open) the output `.msh` in **gmsh**. Show the
field by selecting the **`normE`** view (electric-field magnitude, V/m). The
field on the brain surface is usually what you want — toggle the GM/brain view.

## Notes / limitations of a custom mesh
- **No EEG names** ("Cz", "F3") and **no MNI coordinates** — those need an `m2m`
  folder. Place positions by clicking, as above.
- Default tissue **conductivities apply automatically** because the mesh uses
  SimNIBS standard tissue tags (brain=2, CSF=3, bone=4, scalp=5).
- The **skull is an approximation** (see README) — fine for exploring montages
  and relative comparisons; be cautious about absolute V/m until you add a CT.
