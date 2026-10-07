"""Finer stator meshes to settle convergence: 1.4 and 1.1 mm."""
from paths import *  # repo paths; also switches to the work folder
import sys
import os, subprocess, math
DL = WORK
SOLV = SOLV_DIR
CCX = CCX_EXE
env=dict(os.environ); env['PATH']=SOLV+os.pathsep+env['PATH']; env['OMP_NUM_THREADS']='4'
def sh(c): subprocess.run(c,cwd=DL,env=env,capture_output=True)
prad=lambda th:3452*math.cos(4*th)+4117*math.cos(16*th)
def deck(mesh,faces,out):
    dl=[f'{int(e)}, P{int(f)}, {prad(float(th))*1e-6:.8e}' for e,f,th,_ in
        (ln.split(',') for ln in open(os.path.join(DL,faces)))]
    open(os.path.join(DL,out),'w').write(f'''*Heading
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
*Static
*Dload
'''+'\n'.join(dl)+'\n*Node file\nU\n*End step\n')
def peak(frd):
    mm=0;ib=False
    for ln in open(os.path.join(DL,frd),encoding='latin-1'):
        if ln.startswith(' -4') and 'DISP' in ln: ib=True; continue
        if ib and ln.startswith(' -1'):
            try:
                u=(float(ln[13:25])**2+float(ln[25:37])**2+float(ln[37:49])**2)**.5
                if u>mm:mm=u
            except:pass
        elif ib and ln.startswith(' -3'):break
    return mm*1e6
for s in [1.4,1.1]:
    suf=f'_c{int(s*10)}'; mesh=f'radial_mesh{suf}.inp'; faces=f'radial_faces{suf}.txt'
    sh([sys.executable,os.path.join(SCRIPTS,'radial_mesh.py'),str(s),suf])
    nel=0;ie=False
    for ln in open(os.path.join(DL,mesh)):
        if ln.startswith('*Element'):ie=True;continue
        if ie and ln.startswith('*'):break
        if ie and ln.strip():nel+=1
    job=f'conv_radial{suf}'; deck(mesh,faces,job+'.inp'); sh([CCX,'-i',job])
    print(f'radial size={s}mm elems={nel} peak={peak(job+".frd"):.4f} nm')
