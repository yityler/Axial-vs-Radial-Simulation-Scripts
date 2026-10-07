"""
Order-0 axial forcing of the bore-clamped disc (B1 = baseline: fine mesh with its own
bore-face clamp set FIX, as in disc_harm_fine / conv_disc_c12 -> modes 765 / 876 Hz).
NOTE: bc_study.py 'discB1' clamps r < 10.5 mm instead (modes 788 / 892 Hz); not used here. Three CalculiX runs:
  s0_static : static, uniform 1 kPa           -> static compliance check
  s0_800    : steady-state dynamics at 800 Hz, uniform 12.3 kPa, zeta = 2 %
  s0_sweep  : steady-state dynamics 20-3400 Hz, uniform 1 kPa, zeta = 2 % (unit FRF,
              scaled linearly in speed_sweep_fig.py)
Order 0 = spatially uniform pressure on the loaded disc face.
"""
from paths import *  # repo paths; also switches to the work folder
import math, os, subprocess, sys
DL = WORK
SOLV = SOLV_DIR
CCX = CCX_EXE
env = dict(os.environ); env['PATH'] = SOLV + os.pathsep + env['PATH']; env['OMP_NUM_THREADS'] = '4'

def load_nodes(meshfile):          # same as bc_study.py
    xyz = {}; mode = None
    for ln in open(os.path.join(DL, meshfile)):
        s = ln.strip()
        if s.lower().startswith('*node'): mode = 'n'; continue
        if s.startswith('*') and mode == 'n': break
        if mode == 'n' and s:
            p = s.split(','); xyz[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    return xyz

def load_faces(facefile):          # same as bc_study.py
    out = []
    for ln in open(os.path.join(DL, facefile)):
        e, f, th, r = ln.split(','); out.append((int(e), int(f), float(th)))
    return out

MESH, FACES = 'disc_mesh_fine.inp', 'disc_faces_fine.txt'
ZETA, NMODES = 0.02, 30

xyz = load_nodes(MESH); faces = load_faces(FACES)
rmax = max(math.hypot(x, y) for x, y, z in xyz.values())
rim = [n for n, (x, y, z) in xyz.items() if math.hypot(x, y) > rmax - 0.3]

def nset(name, ids):
    s = f'*Nset, Nset={name}\n'
    for i in range(0, len(ids), 10): s += ', '.join(map(str, ids[i:i + 10])) + ',\n'
    return s

def head(title):
    return (f'*Heading\n{title}\n*Include, Input={MESH}\n*Material, Name=Steel\n*Density\n7.85E-09\n'
            f'*Elastic\n200000, 0.3\n*Solid section, Elset=SOLID, Material=Steel\n'
            + nset('RIM', rim) + '*Boundary\nFIX, 1, 3\n')

def dload(p_pa):
    return '\n'.join(f'{e}, P{f}, {p_pa * 1e-6:.8e}' for e, f, th in faces) + '\n'

def modal():
    return f'*Step\n*Frequency, Solver=Pardiso, Storage=Yes\n{NMODES}\n*End step\n'

decks = {
    's0_static': head('order-0 static 1 kPa') +
        '*Step\n*Static\n*Dload\n' + dload(1000.0) + '*Node file\nU\n*Node print, Nset=RIM\nU\n*End step\n',
    's0_800': head('order-0 SSD 800 Hz 12.3 kPa') + modal() +
        '*Step\n*Steady state dynamics, Solver=Pardiso\n799., 801., 3, 1.\n'
        f'*Modal damping\n1, {NMODES}, {ZETA}\n*Dload, Load case=1\n' + dload(12300.0) +
        '*Node file\nU\n*Node print, Nset=RIM\nU\n*End step\n',
    's0_sweep': head('order-0 SSD sweep unit 1 kPa') + modal() +
        '*Step\n*Steady state dynamics, Solver=Pardiso\n20., 3400., 80, 3.\n'
        f'*Modal damping\n1, {NMODES}, {ZETA}\n*Dload, Load case=1\n' + dload(1000.0) +
        '*Node print, Nset=RIM\nU\n*End step\n',
}
print(f'rim nodes {len(rim)}, loaded faces {len(faces)}, rmax {rmax:.2f} mm', flush=True)
for job in (sys.argv[1:] or decks):
    open(os.path.join(DL, job + '.inp'), 'w').write(decks[job])
    r = subprocess.run([CCX, '-i', job], cwd=DL, env=env, capture_output=True, text=True)
    print(job, 'exit', r.returncode, r.stdout[-300:].replace('\n', ' | '), flush=True)
