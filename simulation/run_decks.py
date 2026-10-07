"""Run the hand-built CalculiX decks in decks/ (they have no generating script).

  disc_modal       coarse disc (disc_mesh.inp), bore-clamped modal       -> fig1b, fig8
  disc_static      fine disc (disc_mesh_fine.inp), static pressure wave  -> fig10, 28.06 um
  disc_frf2        coarse disc, steady-state sweep 50-3000 Hz, rim print -> fig9b
  Analysis-1_wave  radial stator, free-free modal (exported from PrePoMax, self-contained) -> fig1a

Each deck is copied into the work folder and solved there. The disc decks include
the meshes made by disc_mesh.py, so run that first (see README).
Usage: python scripts/run_decks.py [deck ...]     (default: all four)
"""
from paths import *  # repo paths; also switches to the work folder
import os, shutil, subprocess, sys

ALL = ['disc_modal', 'disc_static', 'disc_frf2', 'Analysis-1_wave']
NEEDS = {'disc_modal': 'disc_mesh.inp', 'disc_static': 'disc_mesh_fine.inp', 'disc_frf2': 'disc_mesh.inp'}

env = dict(os.environ); env['PATH'] = SOLV_DIR + os.pathsep + env['PATH']; env.setdefault('OMP_NUM_THREADS', '4')
for job in (sys.argv[1:] or ALL):
    need = NEEDS.get(job)
    if need and not os.path.exists(need):
        sys.exit(f'{job}: {need} not found in {WORK}. Run disc_mesh.py first (see README).')
    shutil.copy2(os.path.join(DECKS, job + '.inp'), job + '.inp')
    r = subprocess.run([CCX_EXE, '-i', job], env=env, capture_output=True, text=True)
    print(f'{job}: exit {r.returncode}' + ('' if r.returncode == 0 else '\n' + r.stdout[-800:]), flush=True)
