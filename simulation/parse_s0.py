"""Parse the order-0 disc runs: peak |U| (complex magnitude, all nodes) from .frd,
and the unit-load rim FRF from s0_sweep.dat."""
from paths import *  # repo paths; also switches to the work folder
import numpy as np

def frd_blocks(path):
    """Yield (freq, name, {node: (u1,u2,u3)}) for every DISP/DISPI block."""
    freq = None; name = None; cur = None
    for ln in open(path, encoding='latin-1'):
        if '100CL' in ln:
            try: freq = float(ln[12:24])
            except ValueError: freq = None
            continue
        if ln.startswith(' -4'):
            name = ln.split()[1]; cur = {} if name in ('DISP', 'DISPI') else None; continue
        if cur is not None:
            if ln.startswith(' -1'):
                cur[int(ln[3:13])] = (float(ln[13:25]), float(ln[25:37]), float(ln[37:49]))
            elif ln.startswith(' -3'):
                yield freq, name, cur; cur = None

def peak_complex(path):
    """{freq: peak over nodes of sqrt(|u1|^2+|u2|^2+|u3|^2)} using real + imaginary parts."""
    re = {}; out = {}
    for f, name, d in frd_blocks(path):
        if name == 'DISP': re[f] = d
        elif name == 'DISPI':
            R = re.pop(f); m = 0.0
            for n, (a, b, c) in R.items():
                i = d.get(n, (0, 0, 0))
                m = max(m, (a*a + b*b + c*c + i[0]**2 + i[1]**2 + i[2]**2) ** 0.5)
            out[f] = m
    return out

def peak_static(path):
    for f, name, d in frd_blocks(path):
        if name == 'DISP':
            return max((a*a + b*b + c*c) ** 0.5 for a, b, c in d.values())

def rim_frf(path):
    """Unit-load FRF: (freq array, peak rim complex magnitude array) from .dat node print."""
    blocks = {}; cur = None
    for ln in open(path, encoding='latin-1'):
        s = ln.strip()
        if s.startswith('displacements') and 'for set RIM' in s:
            f = float(s.split()[-1]); cur = {}; blocks.setdefault(f, []).append(cur); continue
        if cur is not None:
            p = s.split()
            if not p: continue
            try: cur[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            except (ValueError, IndexError): cur = None
    F, U = [], []
    for f in sorted(blocks):
        b = blocks[f]
        if len(b) < 2: continue
        R, I = b[0], b[1]
        U.append(max((sum(x*x for x in R[n]) + sum(x*x for x in I.get(n, (0, 0, 0)))) ** 0.5 for n in R))
        F.append(f)
    return np.array(F), np.array(U)

if __name__ == '__main__':
    st = peak_static('s0_static.frd')
    print(f'static uniform 1 kPa: peak |U| = {st*1e3:.3f} um')
    for f, m in sorted(peak_complex('s0_800.frd').items()):
        print(f'SSD {f:.1f} Hz, 12.3 kPa: peak complex |U| = {m*1e3:.2f} um')
    F, U = rim_frf('s0_sweep.dat')                    # large file; cache what fig12 needs
    np.savez('s0_sweep_frf.npz', F=F, U=U)
    print(f'wrote s0_sweep_frf.npz ({len(F)} frequencies, {F.min():.0f}-{F.max():.0f} Hz)')
