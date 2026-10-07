"""Mesh-convergence study: static peak displacement + first modal freq vs element size,
for disc (bore clamped) and radial stator (outer-rim clamped). Reproducible."""
from paths import *  # repo paths; also switches to the work folder
import sys
import os, subprocess, math, numpy as np
DL = WORK
SOLV = SOLV_DIR
CCX = CCX_EXE
env=dict(os.environ); env['PATH']=SOLV+os.pathsep+env['PATH']; env['OMP_NUM_THREADS']='4'

def sh(cmd): subprocess.run(cmd,cwd=DL,env=env,capture_output=True)

def static_deck(mesh,faces,pfun,out):
    dl=[]
    for ln in open(os.path.join(DL,faces)):
        e,f,th,r=ln.split(','); dl.append(f'{int(e)}, P{int(f)}, {pfun(float(th))*1e-6:.8e}')
    deck=f'''*Heading
conv
*Include, Input={mesh}
*Material, Name=Steel
*Density
7.85E-09
*Elastic
200000, 0.3
*Solid section, Elset=SOLID, Material=Steel
*Boundary
FIX, 1, 3
*Step
*Frequency, Solver=Pardiso
4
*End step
*Step
*Static
*Dload
'''+'\n'.join(dl)+'\n*Node file\nU\n*End step\n'
    open(os.path.join(DL,out),'w').write(deck)

def peak_static(frd):
    mm=0;inblk=False
    for ln in open(os.path.join(DL,frd),encoding='latin-1'):
        if ln.startswith(' -4') and 'DISP' in ln: inblk=True; continue
        if inblk and ln.startswith(' -1'):
            try:
                u1=float(ln[13:25]);u2=float(ln[25:37]);u3=float(ln[37:49]); m=(u1*u1+u2*u2+u3*u3)**.5
                if m>mm:mm=m
            except: pass
        elif inblk and ln.startswith(' -3'): break
    return mm
def first_freq(dat):
    for ln in open(os.path.join(DL,dat),encoding='latin-1'):
        p=ln.split()
        if len(p)>=4 and p[0]=='1':
            try: return float(p[3])
            except: pass

pdisc=lambda th:16400+5347*math.cos(8*th)+5716*math.cos(12*th)
prad =lambda th:3452*math.cos(4*th)+4117*math.cos(16*th)

jobs=[]  # (label, meshscript, size, suffix, pfun, unit)
for s in [3.0,2.0,1.5,1.2]:
    jobs.append(('disc','disc_mesh.py',s,f'_c{int(s*10)}',pdisc,'um'))
for s in [6.0,4.0,2.5,1.8]:
    jobs.append(('radial','radial_mesh.py',s,f'_c{int(s*10)}',prad,'nm'))

rows=[]
for lab,script,size,suf,pfun,unit in jobs:
    mesh=('disc_mesh' if lab=='disc' else 'radial_mesh')+suf+'.inp'
    faces=('disc_faces' if lab=='disc' else 'radial_faces')+suf+'.txt'
    sh([sys.executable,os.path.join(SCRIPTS,script),str(size),suf])
    # element count
    ne=sum(1 for _ in open(os.path.join(DL,mesh))) # rough
    job=f'conv_{lab}{suf}'
    static_deck(mesh,faces,pfun,job+'.inp')
    sh([CCX,'-i',job])
    u=peak_static(job+'.frd'); f1=first_freq(job+'.dat')
    # count real elements
    nel=0; ine=False
    for ln in open(os.path.join(DL,mesh)):
        if ln.startswith('*Element'): ine=True; continue
        if ine and ln.startswith('*'): break
        if ine and ln.strip(): nel+=1
    val=u*1e3 if unit=='um' else u*1e6
    rows.append((lab,size,nel,val,unit,f1))
    print(f'{lab:6s} size={size:4.1f}mm  elems={nel:6d}  peak={val:10.4f} {unit}  f1={f1:.1f} Hz')

np.save('conv_rows.npy',np.array(rows,dtype=object))
print('saved conv_rows.npy')
