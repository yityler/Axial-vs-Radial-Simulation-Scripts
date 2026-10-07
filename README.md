# Axial vs radial flux PMSM vibration (NVH): reproducibility code

Scripts, inputs and figures for a controlled comparison of structural vibration in
an axial-flux and a radial-flux permanent-magnet machine. Pole count, rotating
mass and current density are matched, so the free variable is force direction
(axial "drumhead" loading of a rotor disc vs radial ovalizing of a stator).

The pipeline is CAD (Fusion 360) → magnetic forces (FEMM) → meshing (gmsh) →
modal, static and steady-state dynamic analysis (CalculiX) → figures (matplotlib,
PyVista).

## Repository layout

```
scripts/   all Python scripts (run from the repo root, see below)
data/      magnetic force data from FEMM, cached frequency response, Fusion screenshots
cad/       STEP and STL geometry of the analysed parts and the full radial motor
decks/     hand-built CalculiX input decks that have no generating script
figures/   the final paper figures (PDF + 600-DPI PNG)
work/      created on first run; all meshes, solver files and regenerated figures go here (git-ignored)
```

## Setup

1. Python 3.11 or newer, then
   ```bash
   pip install -r requirements.txt
   ```
2. CalculiX (`ccx`), only needed to rerun the structural analyses. Any recent
   build works. The paper used `ccx_dynamic.exe` from PrePoMax v2.6.0. Point the
   scripts at it with the `CCX` environment variable, or put `ccx` on your PATH.
   ```bash
   export CCX="/path/to/ccx"                 # macOS / Linux / Git Bash
   ```
   ```bat
   set CCX=C:\path\to\PrePoMax\Solver\ccx_dynamic.exe
   ```
3. FEMM 4.2 and `pyfemm`, only needed to rerun the magnetic analyses (Windows).

Every script imports `scripts/paths.py`, which switches to the `work/` folder.
Set `NVH_WORK` to use a different work folder. Run scripts from the repo root as
`python scripts/<name>.py`.

## Quick start: rebuild the figures that need no solver

These read only the committed `data/` and print in seconds.

```bash
python scripts/figs_2_6.py          # figs 3, 4, 5, 6
python scripts/fig2_fix.py          # fig 2 (uses the Fusion screenshots in data/fusion)
python scripts/fig7_convergence.py  # fig 7 (convergence values are written in the script)
python scripts/fig11_airgap.py      # fig 11 (from data/*.csv)
python scripts/fig12_speed_sweep.py # fig 12 (from data/s0_sweep_frf.npz)
```

Output goes to `work/figures/`.

## Full pipeline

Run in this order. Times are for a 4-core laptop.

### 1. Magnetic forces (FEMM, optional, results already in `data/`)

```bash
python scripts/radial_loaded_force_v2.py   # radial machine, loaded, 24 rotor steps -> data/radial_loaded_force.csv
python scripts/axial_slice1.py             # axial machine, single 2D slice (set-up check)
python scripts/axial_multislice.py         # axial machine, multi-slice total axial force (~242 N)
python scripts/axial_sweep.py              # axial machine, rotor sweep -> data/axial_force_sweep.csv
```

### 2. Meshes (gmsh, seconds to about a minute)

```bash
python scripts/disc_mesh.py                # axial rotor disc, 2.0 mm  -> disc_mesh.inp, disc_faces.txt
python scripts/disc_mesh.py 1.2 _fine      # axial rotor disc, 1.2 mm  -> disc_mesh_fine.inp (baseline)
python scripts/radial_mesh.py              # radial stator, 4.0 mm     -> radial_mesh.inp, radial_faces.txt
```

### 3. Structural analyses (CalculiX)

```bash
python scripts/run_decks.py                # hand-built decks: disc_modal, disc_static, disc_frf2, Analysis-1_wave
python scripts/build_disc_harm.py _fine    # disc steady-state deck at 400 Hz
ccx -i disc_harm_fine                      # (run inside work/)
python scripts/build_radial_harm.py        # stator steady-state deck at 200 Hz
ccx -i radial_harm                         # (run inside work/)
python scripts/parse_ssd.py                # peak complex displacement of the two steady-state runs
python scripts/conv_run.py                 # mesh convergence, disc and stator (several minutes)
python scripts/conv_stator.py              # two finer stator meshes (several minutes)
python scripts/bc_study.py                 # support-condition sensitivity (fig 5)
python scripts/speed_sweep_solve.py        # order-0 disc runs for fig 12 (about 8 minutes)
python scripts/parse_s0.py                 # 800 Hz result + writes s0_sweep_frf.npz for fig 12
```

`disc_frf2` and `s0_sweep` write large `.dat` files (about 1 GB and 0.6 GB).

### 4. Figures that read solver output

```bash
python scripts/fig1_modeshapes.py   # needs Analysis-1_wave.frd, disc_modal.frd
python scripts/fig8_validation.py   # needs disc_modal.frd
python scripts/fig9_frf.py          # needs disc_frf2.dat
python scripts/fig10_deflection.py  # needs disc_static.frd
python scripts/render3d.py          # FE-geometry renders, fallback for fig 2
python scripts/render_fullmotor.py  # full radial motor STL render, second fallback for fig 2
```

### 5. Fusion 360 screenshots for fig 2 (optional)

`fshot.py` talks to the MCP server built into Fusion 360 (`127.0.0.1:27184`).
With the motor document open in Fusion:

```bash
python scripts/fshot.py shot data/fusion/_radial3d_fusion.png iso-top-right
```

## Figure index

| Figure | Script | Inputs |
|---|---|---|
| 1 Mode shapes | `fig1_modeshapes.py` | `Analysis-1_wave`, `disc_modal` |
| 2 Geometry | `fig2_fix.py` | `data/fusion/` |
| 3 Force spectra | `figs_2_6.py` | values in script, from the FEMM data |
| 4 Ratio decomposition | `figs_2_6.py` | values in script |
| 5 Support sensitivity | `figs_2_6.py` | values in script, from `bc_study.py` |
| 6 Pipeline | `figs_2_6.py` | none |
| 7 Mesh convergence | `fig7_convergence.py` | values in script, from `conv_run.py`, `conv_stator.py` |
| 8 Analytical validation | `fig8_validation.py` | `disc_modal` |
| 9 Frequency response | `fig9_frf.py` | `disc_frf2` (disc), SDOF model (stator) |
| 10 Static deflection | `fig10_deflection.py` | `disc_static` |
| 11 Air-gap field and stress | `fig11_airgap.py` | `data/*.csv` |
| 12 Speed sweep | `fig12_speed_sweep.py` | `s0_sweep` via `data/s0_sweep_frf.npz` |

Shared plotting style (fonts, sizes, colours, save format) is in `nvh_style.py`.

## Modelling assumptions

* Linear elastic steel, E = 200 GPa, ν = 0.3, ρ = 7850 kg/m³. No housing, winding or magnet mass.
* Units in the CalculiX decks are mm, tonne, s, MPa.
* Modal damping ζ = 2% on all modes.
* Disc: bore face clamped (node set `FIX`). Stator: outer rim clamped.
* Force harmonics are applied with zero relative phase.
* The radial-stator curve in fig 9 is a single-degree-of-freedom model
  (first mode 35,875 Hz, static 1.38 nm), not a finite-element sweep.

## Verified

These were rerun from this repository and compared with the original outputs.

* All meshes and generated decks are byte-identical.
* CalculiX runs reproduce the original eigenfrequencies and peak displacements.
* All twelve figures are pixel-identical.

The FEMM scripts, the long convergence and support-condition runs, and the
Fusion screenshot capture were not rerun.
