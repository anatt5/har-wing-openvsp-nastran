# HAR Wing Sensitivity Studies

Python scripts supporting a thesis on a High-Aspect-Ratio (HAR) wing aircraft: parametric geometry
in **OpenVSP**, aerodynamic/flight-dynamics analysis in **VSPAERO**, and structural (FEA) analysis
in **MSC Nastran**. Each script is standalone and builds its own geometry/mesh from scratch, so the
repo shows the *method* behind each sensitivity study rather than a fixed pipeline you run end to end.

## Structure

```
airfoils/                            Shared airfoil coordinate files (.dat), used by every script
utils/
  airfoil_utils.py                   Shared OpenVSP airfoil-assignment helpers

1_geometry_structure/
  aircraft_geometry.py               Full aircraft geometry: wing + HTP + VTP + fuselage
                                      (loads OpenVSP's stock TransportFuse0.vsp3 template)

2_aero_steady_state/
  panel_density_convergence.py       VSPAERO panel-density (VLM mesh) convergence study
  alpha_sweep.py                     Steady-state alpha sweep -> CL_alpha, CD0, static margin, etc.

3_aero_cg_sweep_dynamics/
  cg_sweep_dynamic_modes.py          CG sweep: trim + Etkin longitudinal/lateral dynamic modes

4_structural_sensitivity/
  mesh_convergence/
    run_mesh_convergence.py          Structural (Nastran) FEA mesh-density convergence study
    save_shared_repaired_mesh.py     Saves one canonical repaired mesh for the other studies to reuse
  taper_sweep/
    run_constant_weight_sweep.py     Skin taper (root/tip thickness) sweep at constant skin mass
  tolerance_sensitivity/
    run_tolerance_sensitivity.py     Node-equivalencing tolerance sensitivity study
  AR_sweep/
    aero/ar_sweep_aero.py            Aspect-ratio sweep -- aerodynamic / flight-dynamics side
    structural/ar_sweep_structural.py  Aspect-ratio sweep -- structural (Nastran) side
```

## Requirements

**Python packages** (`pip install -r requirements.txt`):
- numpy, scipy, matplotlib, pyNastran, plotly

**External software** (not installable via pip):
- **[OpenVSP](https://openvsp.org)**, with its Python API enabled (`import openvsp`) -- required by
  every script.
- **MSC Nastran** (licensed), for the scripts under `4_structural_sensitivity/` -- they call it
  directly via `subprocess.run([NASTRAN_EXE, ...])`. The executable path is hardcoded near the top
  of each script and will need to be updated to match your own installation.

## Running a study

Every script is self-contained: `cd` into its folder and run it directly, e.g.

```bash
cd 4_structural_sensitivity/mesh_convergence
python run_mesh_convergence.py
```

It builds its own geometry/mesh, runs the analysis (VSPAERO or Nastran) for each point in its sweep,
and writes its results (CSV, figures, log) into that same folder.

Most scripts are fully independent and can be run in any order:
`1_geometry_structure/aircraft_geometry.py`, both scripts in `2_aero_steady_state/`,
`3_aero_cg_sweep_dynamics/cg_sweep_dynamic_modes.py`, and both scripts in `4_structural_sensitivity/AR_sweep/`.

The four scripts inside `4_structural_sensitivity/` that share a mesh, however, **must be run in
this order**, since each one reads a file produced by the previous step:

1. `mesh_convergence/run_mesh_convergence.py` -- builds the raw mesh for every BaseLen, including
   `BaseLen_0.50/AC_wing.bdf` (the thesis-wide reference density).
2. `mesh_convergence/save_shared_repaired_mesh.py` -- reads that `BaseLen_0.50/AC_wing.bdf` and
   writes the repaired snapshot `AC_wing_repaired_half_wing.bdf`.
3. `taper_sweep/run_constant_weight_sweep.py` -- reads the repaired snapshot from step 2.
4. `tolerance_sensitivity/run_tolerance_sensitivity.py` -- reads the *raw* `BaseLen_0.50/AC_wing.bdf`
   from step 1 directly (it repairs it itself, once per candidate tolerance) -- does **not** need
   step 2.

Running 2, 3, or 4 before the file they depend on exists will fail with a missing-file error.

## Shared code

The airfoil-assignment logic used by every geometry-building script lives once in
`utils/airfoil_utils.py` instead of being duplicated per script.
