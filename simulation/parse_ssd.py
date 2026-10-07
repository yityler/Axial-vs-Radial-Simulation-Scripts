"""Correct steady-state-dynamics .frd parser: combine DISP (real) + DISPI (imag)
into complex displacement amplitude per node, report peak at a target frequency."""
from paths import *  # repo paths; also switches to the work folder
import math, sys

def peak_ssd(fn, ftarget, ftol=0.6):
    # collect, per (freq, kind), node -> (u1,u2,u3)
    freq=None; kind=None; real={}; imag={}
    def store(nid,vals):
        (real if kind=='DISP' else imag).setdefault(freq,{})[nid]=vals
    for ln in open(fn,encoding='latin-1'):
        if '100CL' in ln:
            try: freq=float(ln.split()[2])
            except: freq=None
            kind=None
        elif ln.startswith(' -4'):
            kind='DISPI' if 'DISPI' in ln else ('DISP' if 'DISP' in ln else None)
        elif ln.startswith(' -1') and freq is not None and kind in ('DISP','DISPI') and len(ln)>48:
            try:
                nid=int(ln[3:13]); u1=float(ln[13:25]); u2=float(ln[25:37]); u3=float(ln[37:49])
                store(nid,(u1,u2,u3))
            except: pass
    # find the frequency closest to target present in both
    freqs=[f for f in real if any(abs(f-ftarget)<ftol for _ in [0])]
    fbest=min((f for f in real if abs(f-ftarget)<ftol), key=lambda x:abs(x-ftarget), default=None)
    if fbest is None: return None
    R=real[fbest]; I=imag.get(fbest,{})
    mm=0; mmU3=0
    for nid,(r1,r2,r3) in R.items():
        i1,i2,i3 = I.get(nid,(0,0,0))
        mag=math.sqrt(r1*r1+i1*i1 + r2*r2+i2*i2 + r3*r3+i3*i3)
        u3amp=math.sqrt(r3*r3+i3*i3)
        if mag>mm: mm=mag
        if u3amp>mmU3: mmU3=u3amp
    return fbest, mm, mmU3

for fn,ft,lbl in [('disc_harm_fine.frd',400.0,'DISC 1.2mm @400Hz'),
                  ('disc_harm.frd',400.0,'DISC 2.0mm @400Hz'),
                  ('radial_harm.frd',200.0,'STATOR 4mm @200Hz')]:
    try:
        f,mm,u3=peak_ssd(fn,ft)
        unit='um' if mm>1e-2 else 'nm'
        val=mm*1e3 if unit=='um' else mm*1e6
        print(f'{lbl:22s}: peak|U|(complex) = {mm*1e6:.3f} nm = {mm*1e3:.4f} um   (max|U3|={u3*1e6:.1f} nm)')
    except Exception as e:
        print(f'{lbl}: {e}')
