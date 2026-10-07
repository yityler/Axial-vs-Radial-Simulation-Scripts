"""
Study A (compliance/amplification decomposition) + Study B (BC sensitivity).
Same gmsh+ccx pipeline. u_dyn = u_static * SDOF amplification 1/(1-(f/f1)^2),
because modal steady-state-dynamics is unreliable below the first mode (verified).
All values PROVISIONAL.
"""
from paths import *  # repo paths; also switches to the work folder
import math, os, subprocess, sys

DL = WORK
SOLV = SOLV_DIR
CCX = CCX_EXE
env = dict(os.environ); env['PATH'] = SOLV + os.pathsep + env['PATH']; env['OMP_NUM_THREADS']='4'

def load_nodes(meshfile):
    xyz={}; mode=None
    for ln in open(os.path.join(DL,meshfile)):
        s=ln.strip()
        if s.lower().startswith('*node'): mode='n'; continue
        if s.startswith('*') and mode=='n': break
        if mode=='n' and s:
            p=s.split(','); xyz[int(p[0])]=(float(p[1]),float(p[2]),float(p[3]))
    return xyz

def load_faces(facefile):
    out=[]
    for ln in open(os.path.join(DL,facefile)):
        e,f,th,r=ln.split(','); out.append((int(e),int(f),float(th)))
    return out

def pfun_disc(th):  return 16400.0 + 5347.0*math.cos(8*th) + 5716.0*math.cos(12*th)
def pfun_rad(th):   return 3452.0*math.cos(4*th) + 4117.0*math.cos(16*th)

def run_variant(name, meshfile, facefile, pfun, fixrule, forcing):
    xyz=load_nodes(meshfile); faces=load_faces(facefile)
    fix=[nid for nid,(x,y,z) in xyz.items() if fixrule(math.hypot(x,y), math.atan2(y,x), z)]
    dl=[f'{e}, P{f}, {pfun(th)*1e-6:.8e}' for e,f,th in faces]
    deck=f'''*Heading
{name}
*Include, Input={meshfile}
*Material, Name=Steel
*Density
7.85E-09
*Elastic
200000, 0.3
*Solid section, Elset=SOLID, Material=Steel
*Nset, Nset=MYFIX
'''
    for i in range(0,len(fix),10): deck+=', '.join(map(str,fix[i:i+10]))+',\n'
    deck+='''*Boundary
MYFIX, 1, 3
*Step
*Frequency, Solver=Pardiso
8
*End step
*Step
*Static
*Dload
'''
    deck+='\n'.join(dl)+'\n*Node file\nU\n*End step\n'
    job=f'bc_{name}'
    open(os.path.join(DL,job+'.inp'),'w').write(deck)
    subprocess.run([CCX,'-i',job],cwd=DL,env=env,capture_output=True)
    # first eigenfrequency from .dat
    f1=None
    for ln in open(os.path.join(DL,job+'.dat'),encoding='latin-1'):
        p=ln.split()
        if len(p)>=4 and p[0]=='1':
            try: f1=float(p[3]); break
            except: pass
    # peak |U| static from .frd
    mm=0; inblk=False
    for ln in open(os.path.join(DL,job+'.frd'),encoding='latin-1'):
        if ln.startswith(' -4') and 'DISP' in ln: inblk=True; continue
        if inblk and ln.startswith(' -1'):
            try:
                u1=float(ln[13:25]);u2=float(ln[25:37]);u3=float(ln[37:49])
                m=(u1*u1+u2*u2+u3*u3)**.5
                if m>mm: mm=m
            except: pass
        elif inblk and ln.startswith(' -3'): break
    Ppk=max(abs(pfun(2*math.pi*k/720)) for k in range(720))
    amp=1.0/abs(1-(forcing/f1)**2) if f1 else float('nan')
    return dict(name=name,nfix=len(fix),f1=f1,forcing=forcing,fratio=forcing/f1,
                u_static=mm*1e6,  # nm
                Ppk=Ppk, amp=amp, u_dyn=mm*1e6*amp)

PAD=lambda th,centers,half: any(abs((th-c+math.pi)%(2*math.pi)-math.pi)<half for c in centers)
c4=[0,math.pi/2,math.pi,3*math.pi/2]

variants=[
 # disc: fine mesh, forcing 400 Hz
 ('discB1','disc_mesh_fine.inp','disc_faces_fine.txt',pfun_disc, lambda r,th,z: r<10.5, 400.0),
 ('discB2','disc_mesh_fine.inp','disc_faces_fine.txt',pfun_disc, lambda r,th,z: r<20.0, 400.0),
 ('discB3','disc_mesh_fine.inp','disc_faces_fine.txt',pfun_disc, lambda r,th,z: r<20.0 or r>62.0, 400.0),
 # stator: 4mm mesh, forcing 200 Hz
 ('statS1','radial_mesh.inp','radial_faces.txt',pfun_rad, lambda r,th,z: abs(r-50.0)<0.8, 200.0),
 ('statS2','radial_mesh.inp','radial_faces.txt',pfun_rad, lambda r,th,z: abs(r-50.0)<0.8 and PAD(th,c4,math.radians(15)), 200.0),
]

res=[run_variant(*v) for v in variants]
print(f'{"variant":8s} {"nfix":>6s} {"f1(Hz)":>10s} {"f/f1":>8s} {"u_stat(nm)":>12s} {"Ppk(N/m2)":>10s} {"amp":>7s} {"u_dyn(nm)":>12s}')
for r in res:
    print(f'{r["name"]:8s} {r["nfix"]:6d} {r["f1"]:10.1f} {r["fratio"]:8.4f} {r["u_static"]:12.3f} {r["Ppk"]:10.0f} {r["amp"]:7.3f} {r["u_dyn"]:12.3f}')

d={r['name']:r for r in res}
print('\n--- STUDY A decomposition (baseline discB1 vs statS1) ---')
D,S=d['discB1'],d['statS1']
force=D['Ppk']/S['Ppk']
compl=(D['u_static']/D['Ppk'])/(S['u_static']/S['Ppk'])
ampf=D['amp']/S['amp']
print(f'  force factor      P_disc/P_stat        = {force:.2f}x')
print(f'  compliance factor (u/P)_disc/(u/P)_stat = {compl:.1f}x')
print(f'  amplification     amp_disc/amp_stat     = {ampf:.2f}x')
print(f'  product                                 = {force*compl*ampf:.0f}x')
print(f'  direct u_dyn ratio disc/stat            = {D["u_dyn"]/S["u_dyn"]:.0f}x')

print('\n--- STUDY B matched-BC ratios (u_dyn disc/stator) ---')
for db in ('discB2','discB3'):
    print(f'  {db} / statS2 = {d[db]["u_dyn"]/d["statS2"]["u_dyn"]:.0f}x   (baseline discB1/statS1 = {D["u_dyn"]/S["u_dyn"]:.0f}x)')
