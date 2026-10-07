# Simulation scripts

Everything needed to rerun the analyses, with no figure code. Run from the repo
root. Output goes to `work/` (see the main README for setup, `CCX` and `NVH_WORK`).
The figure scripts that read these results are in `../plotting/`.

## 1. Magnetic forces (FEMM 4.2 + pyfemm, Windows)

Optional. Their results are already committed in `data/`.

| Script | What it does | Output |
|---|---|---|
| `radial_loaded_force_v2.py` | Radial machine at rated load (6 A/mm², 3000 r/min), 24 rotor steps, air-gap Maxwell stress | `data/radial_loaded_force.csv` |
| `axial_slice1.py` | Axial machine, single 2D slice, set-up check | printed |
| `axial_multislice.py` | Axial machine, several radial slices, total axial force (~242 N) | printed |
| `axial_sweep.py` | Axial machine, rotor sweep at the mean-radius slice | `data/axial_force_sweep.csv` |

## 2. Meshes (gmsh)

| Command | Output |
|---|---|
| `python simulation/disc_mesh.py` | `disc_mesh.inp`, `disc_faces.txt` (2.0 mm) |
| `python simulation/disc_mesh.py 1.2 _fine` | `disc_mesh_fine.inp`, `disc_faces_fine.txt` (1.2 mm, baseline) |
| `python simulation/radial_mesh.py` | `radial_mesh.inp`, `radial_faces.txt` (4.0 mm) |

## 3. Structural analyses (CalculiX)

Run in this order.

| Command | What it does |
|---|---|
| `python simulation/run_decks.py` | Hand-built decks from `decks/`: disc modal, disc static, disc frequency sweep, stator free-free modal |
| `python simulation/build_disc_harm.py _fine` then `ccx -i disc_harm_fine` | Disc steady-state response at 400 Hz |
| `python simulation/build_radial_harm.py` then `ccx -i radial_harm` | Stator steady-state response at 200 Hz |
| `python simulation/parse_ssd.py` | Peak complex displacement of the two runs above |
| `python simulation/conv_run.py` | Mesh convergence, disc and stator |
| `python simulation/conv_stator.py` | Two finer stator meshes |
| `python simulation/bc_study.py` | Support-condition sensitivity |
| `python simulation/speed_sweep_solve.py` | Order-0 disc runs: static, 800 Hz, 20 to 3400 Hz sweep |
| `python simulation/parse_s0.py` | 800 Hz peak, and writes `s0_sweep_frf.npz` for the speed-sweep figure |

Run the two `ccx -i` commands inside `work/`. Everything else can be run from
the repo root.

## Helpers

| File | Role |
|---|---|
| `paths.py` | Repo paths, solver location, switches to `work/` |
| `parse_ssd.py` | Steady-state `.frd` reader (real + imaginary parts) |
| `parse_s0.py` | Readers for the order-0 runs (`.frd` peak and rim `.dat` sweep) |
