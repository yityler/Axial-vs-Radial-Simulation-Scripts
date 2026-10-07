"""Repo-relative paths, shared by every script.

Importing this module also changes the working directory to the work folder,
so the scripts' relative file names (meshes, solver decks, results, figures/)
all land in one git-ignored place.

  NVH_WORK  work folder (default: <repo>/work)
  CCX       CalculiX executable (default: 'ccx' on PATH). Example on Windows:
            set CCX=C:\\path\\to\\PrePoMax\\Solver\\ccx_dynamic.exe
"""
import os

SCRIPTS = os.path.dirname(os.path.abspath(__file__))   # this folder (simulation or plotting)
REPO = os.path.dirname(SCRIPTS)
DATA = os.path.join(REPO, 'data')
CAD = os.path.join(REPO, 'cad')
DECKS = os.path.join(REPO, 'decks')
WORK = os.path.abspath(os.environ.get('NVH_WORK', os.path.join(REPO, 'work')))
CCX_EXE = os.environ.get('CCX', 'ccx')
SOLV_DIR = os.path.dirname(CCX_EXE) or os.curdir

os.makedirs(os.path.join(WORK, 'figures'), exist_ok=True)
os.chdir(WORK)
