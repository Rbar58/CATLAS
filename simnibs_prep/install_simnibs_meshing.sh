#!/usr/bin/env bash
# Install ONLY the SimNIBS meshing capability (create_mesh / meshmesh) on a
# headless Linux x86_64 box with system Python 3.11 -- no conda, no GUI, no GPU.
#
# SimNIBS is normally installed via its own self-contained installer and is not
# on PyPI. This grabs the official cp311 Linux wheel directly and installs just
# enough to run the CGAL mesher, stubbing the human-only / GPU-only pieces
# (charm's torch & brainnet, the FEM solver's mumps) that meshing never calls.
#
# This is what made cloud-side meshing of the cat model possible. The FEM
# *simulation* still needs the full SimNIBS install on your own machine.
set -euo pipefail
V=4.6.0
BASE="https://github.com/simnibs/simnibs/releases/download/v${V}"

# 1. SimNIBS wheel (no deps -- it pulls non-PyPI packages we handle ourselves)
pip install --no-deps "${BASE}/simnibs-${V}-cp311-cp311-linux_x86_64.whl"

# 2. ordinary PyPI deps (+ antspyx for the accurate Stage-1 registration, and
#    SimpleITK for the fallback registration)
pip install h5py jsonschema pillow requests numba gmsh scipy nibabel scikit-learn SimpleITK antspyx

# 3. SimNIBS-channel wheels (compiled, cp311 manylinux)
pip install --no-deps \
  "https://github.com/simnibs/cortech/releases/download/v0.1/cortech-0.1-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl" \
  "https://github.com/simnibs/fmm3dpy/releases/download/v1.0.4/fmm3dpy-1.0.4-cp311-cp311-manylinux_2_28_x86_64.whl" \
  "https://github.com/simnibs/petsc4py/releases/download/v3.22.2/petsc4py-3.22.2-cp311-cp311-manylinux_2_28_x86_64.whl"

# 4. runtime native libs (petsc4py needs MKL .so.2; CGAL ext needs libtbb)
pip install "mkl==2024.2.2" tbb pygpc

# 5. stub the modules only the (unused) FEM solver / surface-NN code imports
SP=$(python3 -c "import site; print(site.getsitepackages()[0])")
cat > "$SP/mumps.py" <<'PY'
class Context:
    def __init__(self,*a,**k): raise RuntimeError("mumps stub: FEM solver not installed (meshing needs none).")
def __getattr__(n): raise RuntimeError(f"mumps stub: '{n}' unavailable.")
PY
cat > "$SP/torch.py" <<'PY'
class device:
    def __init__(self,*a,**k): pass
class Tensor: pass
def tensor(*a,**k): raise RuntimeError("torch stub: only needed for surface-NN, not meshing.")
def __getattr__(n): raise RuntimeError(f"torch stub: '{n}' unavailable.")
PY
mkdir -p "$SP/brainnet"
cat > "$SP/brainnet/__init__.py" <<'PY'
from . import helpers, datasets
class Surface: pass
class DeepSurferTopology: pass
def __getattr__(n): return type(n,(object,),{})
PY
cat > "$SP/brainnet/helpers.py"  <<'PY'
def __getattr__(n): return type(n,(object,),{})
PY
cat > "$SP/brainnet/datasets.py" <<'PY'
import numpy as _np
class _M:
    def numpy(self): return _np.eye(4, dtype=_np.float64)
MNI305_to_MNI152 = _M()
def __getattr__(n): return type(n,(object,),{})
PY

# 6. verify (MKL/TBB live in this prefix's lib dir)
export LD_LIBRARY_PATH="$(python3 -c 'import sys,os;print(os.path.join(sys.prefix,"lib"))'):${LD_LIBRARY_PATH:-}"
python3 -c "from simnibs.mesh_tools.meshing import create_mesh; import simnibs; print('OK SimNIBS', simnibs.__version__)"
echo "Meshing ready. Remember to export LD_LIBRARY_PATH (see step 6) before meshing."
