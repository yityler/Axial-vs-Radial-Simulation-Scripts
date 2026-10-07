"""
RADIAL STATOR - mesh + CalculiX fragments (gmsh, same pipeline as disc_mesh.py)
===============================================================================
Rebuilds the radial stator through gmsh + CalculiX so it uses IDENTICAL tooling
to the axial disc, for an airtight same-tooling, both-converged comparison.

Writes:
  radial_mesh{SUF}.inp  : *Node, *Element(C3D10), *Nset FIX (outer rim r~50),
                          *Elset SOLID
  radial_faces{SUF}.txt : loaded tooth-tip faces  elem,facenum,theta_rad,r_mm
                          (bore-facing, radial-normal faces, all axial z)

Usage: python radial_mesh.py [elsize_mm] [suffix]
Units: mm / tonne / s / MPa (matches radial + axial structural decks).
"""
from paths import *  # repo paths; also switches to the work folder
import os
import gmsh, math, sys

STEP = os.path.join(CAD, 'Radial_Stator_modal.stp')
ELSIZE = float(sys.argv[1]) if len(sys.argv)>1 else 4.0
SUF    = sys.argv[2] if len(sys.argv)>2 else ''
R_OUT  = 50.0          # outer rim radius (mm) -> fixed
RTOL   = 0.8           # tol for outer-rim node pick (mm)

gmsh.initialize(); gmsh.option.setNumber('General.Terminal',0)
gmsh.open(STEP)
gmsh.option.setNumber('Mesh.MeshSizeMin', ELSIZE*0.6)
gmsh.option.setNumber('Mesh.MeshSizeMax', ELSIZE)
gmsh.option.setNumber('Mesh.ElementOrder', 2)
gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 0)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
gmsh.model.mesh.generate(3)

# nodes (renumber 1..N)
ntags,ncoords,_ = gmsh.model.mesh.getNodes(); ncoords=ncoords.reshape(-1,3)
tag2new={}; nodes=[]
for i,t in enumerate(ntags):
    tag2new[int(t)]=i+1; x,y,z=ncoords[i]; nodes.append((i+1,x,y,z))
newid2xyz={n[0]:(n[1],n[2],n[3]) for n in nodes}

# C3D10 elements with authoritative gmsh->Abaqus permutation
_,_,_,nn10,lc,_ = gmsh.model.mesh.getElementProperties(11)
P=[(lc[3*i],lc[3*i+1],lc[3*i+2]) for i in range(nn10)]
cpar={0:(0,0,0),1:(1,0,0),2:(0,1,0),3:(0,0,1)}
def classify(pt,tol=1e-6):
    for i,cp in cpar.items():
        if all(abs(pt[k]-cp[k])<tol for k in range(3)): return ('c',i)
    for i in range(4):
        for j in range(i+1,4):
            mp=tuple((cpar[i][k]+cpar[j][k])/2 for k in range(3))
            if all(abs(pt[k]-mp[k])<tol for k in range(3)): return ('m',(i,j))
    raise RuntimeError('unmatched')
loc=[classify(p) for p in P]
ABQ=[('c',0),('c',1),('c',2),('c',3),('m',(0,1)),('m',(1,2)),('m',(0,2)),
     ('m',(0,3)),('m',(1,3)),('m',(2,3))]
def fe(tok):
    if tok[0]=='c': return loc.index(tok)
    a,b=tok[1]
    for idx,l in enumerate(loc):
        if l[0]=='m' and set(l[1])==set((a,b)): return idx
    raise RuntimeError('edge')
PERM=[fe(t) for t in ABQ]

etypes,etags,enodes=gmsh.model.mesh.getElements(3)
elems=[]; eid=0
for et,tags,conn in zip(etypes,etags,enodes):
    if et!=11: continue
    conn=conn.reshape(-1,10)
    for row in conn:
        eid+=1
        nn=[tag2new[int(t)] for t in row]
        elems.append((eid,[nn[p] for p in PERM]))

def svol6(nn):
    a=newid2xyz[nn[0]];b=newid2xyz[nn[1]];c=newid2xyz[nn[2]];d=newid2xyz[nn[3]]
    ab=[b[i]-a[i] for i in range(3)];ac=[c[i]-a[i] for i in range(3)];ad=[d[i]-a[i] for i in range(3)]
    cx=ac[1]*ad[2]-ac[2]*ad[1];cy=ac[2]*ad[0]-ac[0]*ad[2];cz=ac[0]*ad[1]-ac[1]*ad[0]
    return ab[0]*cx+ab[1]*cy+ab[2]*cz
SWAP=[0,1,3,2,4,8,7,6,5,9]; nflip=0
elems=[(e,(nn if svol6(nn)>0 else [nn[i] for i in SWAP])) for e,nn in elems]

# FIX: outer rim nodes (r ~ R_OUT)
fix=[nid for nid,x,y,z in nodes if abs(math.hypot(x,y)-R_OUT)<RTOL]

# loaded tooth-tip faces: boundary tris with radial normal near the bore radius
FACE_CORNERS={frozenset((0,1,2)):1,frozenset((0,1,3)):2,frozenset((1,2,3)):3,frozenset((0,2,3)):4}
face_lut={}
for e,nn in elems:
    c=nn[:4]
    for locset,fnum in FACE_CORNERS.items():
        face_lut[frozenset(c[i] for i in locset)]=(e,fnum)

bt,btags,bconn=gmsh.model.mesh.getElements(2)
cand=[]   # (elem,face,theta,r, |n.rhat|)
for et,tags,conn in zip(bt,btags,bconn):
    if et!=9: continue
    conn=conn.reshape(-1,6)
    for row in conn:
        cor=[tag2new[int(t)] for t in row[:3]]
        p0,p1,p2=[newid2xyz[c] for c in cor]
        # normal
        u=[p1[i]-p0[i] for i in range(3)]; v=[p2[i]-p0[i] for i in range(3)]
        nx=u[1]*v[2]-u[2]*v[1]; ny=u[2]*v[0]-u[0]*v[2]; nz=u[0]*v[1]-u[1]*v[0]
        nl=math.sqrt(nx*nx+ny*ny+nz*nz) or 1.0
        nx,ny,nz=nx/nl,ny/nl,nz/nl
        cx=(p0[0]+p1[0]+p2[0])/3; cy=(p0[1]+p1[1]+p2[1])/3
        r=math.hypot(cx,cy) or 1e-9
        rhatx,rhaty=cx/r,cy/r
        ndotr=abs(nx*rhatx+ny*rhaty)
        if abs(nz)<0.3 and ndotr>0.7:   # axis-parallel, radial-facing
            hit=face_lut.get(frozenset(cor))
            if hit: cand.append((hit[0],hit[1],math.atan2(cy,cx),r,ndotr))
# bore band = radial faces at the smallest radii (tooth tips)
rmin=min(c[3] for c in cand)
loaded=[c for c in cand if c[3] < rmin+2.0]

with open(rf'radial_mesh{SUF}.inp','w') as f:
    f.write('*Node\n')
    for nid,x,y,z in nodes: f.write(f'{nid}, {x:.6f}, {y:.6f}, {z:.6f}\n')
    f.write('*Element, Type=C3D10, Elset=SOLID\n')
    for e,nn in elems: f.write(f'{e}, '+', '.join(map(str,nn))+'\n')
    f.write('*Nset, Nset=FIX\n')
    for i in range(0,len(fix),10): f.write(', '.join(map(str,fix[i:i+10]))+',\n')
with open(rf'radial_faces{SUF}.txt','w') as f:
    for e,fn,th,r,_ in loaded: f.write(f'{e},{fn},{th:.6f},{r:.4f}\n')

print(f'nodes={len(nodes)}  C3D10={len(elems)}')
print(f'FIX(outer rim) nodes={len(fix)}   tooth-tip loaded faces={len(loaded)}')
print(f'bore tooth-tip radius ~ {rmin:.2f} mm (band to {rmin+2.0:.2f})')
gmsh.finalize()
