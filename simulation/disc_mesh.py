"""
AXIAL ROTOR DISC - mesh + CalculiX deck generator   BOX 4
=========================================================
Meshes Axial_RotorDisc_modal.stp (shaft-bore-clamped annular plate) with
second-order tets (C3D10), matching the radial structural pipeline. Writes:
  disc_mesh.inp   : *Node, *Element(C3D10), *Nset FIX (bore), *Elset SOLID
  disc_faces.txt  : loaded element faces on the flat disc face, active band,
                    as  elem,facenum,theta_rad,r_mm  (for the pressure wave)

Units land in mm / tonne / s / MPa (same consistent set as the radial deck),
so density 7.85e-9, E=200000 MPa give Hz and mm.

FIRST-PASS SIMPLIFICATIONS (flagged):
  - bare steel disc, no magnet/winding mass
  - element size ELSIZE (coarse-ish first pass), convergence still TODO
  - fixed bore (shaft clamp); free outer edge
  - load applied on ONE flat face over the active band r in [R_ACT_IN, R_out]
"""
from paths import *  # repo paths; also switches to the work folder
import os
import gmsh, math, sys

STEP = os.path.join(CAD, 'Axial_RotorDisc_modal.stp')
ELSIZE = float(sys.argv[1]) if len(sys.argv)>1 else 2.0   # target element size (mm)
SUF = sys.argv[2] if len(sys.argv)>2 else ''              # output filename suffix
R_ACT_IN = 21.0        # inner radius of active (magnet) band (mm) - from magnetic model
BORE_R = 9.55          # shaft bore radius (mm)
LOAD_Z_TOP = True      # apply on the top flat face (z max) of the back iron
RTOL = 0.4             # radial tolerance for picking bore nodes (mm)

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.open(STEP)
gmsh.option.setNumber('Mesh.MeshSizeMin', ELSIZE*0.6)
gmsh.option.setNumber('Mesh.MeshSizeMax', ELSIZE)
gmsh.option.setNumber('Mesh.ElementOrder', 2)
gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 0)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
gmsh.model.mesh.generate(3)

# ---- nodes (renumber to 1..N) ----
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
ncoords = ncoords.reshape(-1,3)
tag2new = {}
nodes = []   # (newid, x,y,z)
for i,t in enumerate(ntags):
    nid = i+1
    tag2new[int(t)] = nid
    x,y,z = ncoords[i]
    nodes.append((nid,x,y,z))
newid2xyz = {n[0]:(n[1],n[2],n[3]) for n in nodes}

# ---- C3D10 elements (gmsh tet10 order == Abaqus C3D10 order) ----
etypes, etags, enodes = gmsh.model.mesh.getElements(3)
elems = []   # (eid, [10 node ids])
eid = 0
for et, tags, conn in zip(etypes, etags, enodes):
    if et != 11:      # 11 = 10-node tetrahedron
        continue
    conn = conn.reshape(-1,10)
    for row in conn:
        eid += 1
        elems.append((eid, [tag2new[int(t)] for t in row]))
if not elems:
    print("ERROR: no C3D10 elements generated"); sys.exit(1)

# ---- derive authoritative gmsh->Abaqus C3D10 node permutation ----
# gmsh gives parametric coords of each local node; match midsides to edge midpoints.
_,_,_,nn10,lc,_ = gmsh.model.mesh.getElementProperties(11)
P = [(lc[3*i],lc[3*i+1],lc[3*i+2]) for i in range(nn10)]     # 10 param points
corner_par = {0:(0.,0.,0.),1:(1.,0.,0.),2:(0.,1.,0.),3:(0.,0.,1.)}
def which(pt,tol=1e-6):
    # identify param point: a corner (return ('c',i)) or edge midpoint (('m',(i,j)))
    for i,cp in corner_par.items():
        if all(abs(pt[k]-cp[k])<tol for k in range(3)): return ('c',i)
    for i in range(4):
        for j in range(i+1,4):
            mp=tuple((corner_par[i][k]+corner_par[j][k])/2 for k in range(3))
            if all(abs(pt[k]-mp[k])<tol for k in range(3)): return ('m',(i,j))
    raise RuntimeError(f'unmatched param pt {pt}')
loc = [which(p) for p in P]
def find(tok):
    return loc.index(tok)
# Abaqus C3D10 slots: corners 1-4, then m(1,2),m(2,3),m(3,1),m(1,4),m(2,4),m(3,4)  (0-based edges below)
ABQ = [('c',0),('c',1),('c',2),('c',3),
       ('m',(0,1)),('m',(1,2)),('m',(0,2)),('m',(0,3)),('m',(1,3)),('m',(2,3))]
def findedge(tok):
    if tok[0]=='c': return loc.index(tok)
    a,b=tok[1]
    for idx,l in enumerate(loc):
        if l[0]=='m' and set(l[1])==set((a,b)): return idx
    raise RuntimeError('edge not found')
PERM = [findedge(t) for t in ABQ]
print(f'gmsh->Abaqus C3D10 perm = {PERM}')
elems = [(eid_, [nn[p] for p in PERM]) for eid_,nn in elems]

# orientation check (positive signed volume required by CalculiX)
def svol6(nn):
    a=newid2xyz[nn[0]]; b=newid2xyz[nn[1]]; c=newid2xyz[nn[2]]; d=newid2xyz[nn[3]]
    ab=(b[0]-a[0],b[1]-a[1],b[2]-a[2]); ac=(c[0]-a[0],c[1]-a[1],c[2]-a[2]); ad=(d[0]-a[0],d[1]-a[1],d[2]-a[2])
    cx=ac[1]*ad[2]-ac[2]*ad[1]; cy=ac[2]*ad[0]-ac[0]*ad[2]; cz=ac[0]*ad[1]-ac[1]*ad[0]
    return ab[0]*cx+ab[1]*cy+ab[2]*cz
SWAP=[0,1,3,2,4,8,7,6,5,9]
nflip=0; fixed=[]
for eid_,nn in elems:
    if svol6(nn)<0: nn=[nn[i] for i in SWAP]; nflip+=1
    fixed.append((eid_,nn))
elems=fixed
print(f'orientation: flipped {nflip}/{len(elems)} tets')

# ---- bore node set (fixed): nodes with radius ~ BORE_R ----
fix_nodes = []
for nid,x,y,z in nodes:
    r = math.hypot(x,y)
    if abs(r-BORE_R) < RTOL:
        fix_nodes.append(nid)

# ---- loaded faces: boundary triangles on the chosen flat face, active band ----
# CalculiX C3D10 face corner-node sets -> face id
FACE_CORNERS = {frozenset((0,1,2)):1, frozenset((0,1,3)):2,
                frozenset((1,2,3)):3, frozenset((0,2,3)):4}   # 0-based local indices
# build lookup: frozenset(global corner ids) -> (eid, facenum)
face_lut = {}
for eid_,nn in elems:
    c = nn[:4]
    for locset,fnum in FACE_CORNERS.items():
        key = frozenset(c[i] for i in locset)
        face_lut[key] = (eid_,fnum)

# z of the flat face to load
zc = [z for _,_,_,z in nodes]
zface = max(zc) if LOAD_Z_TOP else min(zc)
ztol = 0.25

# 2D boundary triangles from gmsh (type 9 = 6-node tri); use their 3 corners
bt, btags, bconn = gmsh.model.mesh.getElements(2)
loaded = []
seen = set()
for et, tags, conn in zip(bt, btags, bconn):
    if et != 9:      # 9 = 6-node triangle
        continue
    conn = conn.reshape(-1,6)
    for row in conn:
        corners = [tag2new[int(t)] for t in row[:3]]
        xyz = [newid2xyz[c] for c in corners]
        zc3 = [p[2] for p in xyz]
        if not all(abs(zz-zface) < ztol for zz in zc3):
            continue
        cx = sum(p[0] for p in xyz)/3; cy = sum(p[1] for p in xyz)/3
        r = math.hypot(cx,cy)
        if r < R_ACT_IN:
            continue
        key = frozenset(corners)
        hit = face_lut.get(key)
        if hit is None:
            continue
        if hit in seen:
            continue
        seen.add(hit)
        theta = math.atan2(cy,cx)
        loaded.append((hit[0],hit[1],theta,r))

# ---- write mesh inp ----
with open(rf'disc_mesh{SUF}.inp','w') as f:
    f.write('*Node\n')
    for nid,x,y,z in nodes:
        f.write(f'{nid}, {x:.6f}, {y:.6f}, {z:.6f}\n')
    f.write('*Element, Type=C3D10, Elset=SOLID\n')
    for eid_,nn in elems:
        f.write(f'{eid_}, ' + ', '.join(str(n) for n in nn) + '\n')
    f.write('*Nset, Nset=FIX\n')
    for i in range(0,len(fix_nodes),10):
        f.write(', '.join(str(n) for n in fix_nodes[i:i+10]) + ',\n')

with open(rf'disc_faces{SUF}.txt','w') as f:
    for eid_,fnum,theta,r in loaded:
        f.write(f'{eid_},{fnum},{theta:.6f},{r:.4f}\n')

print(f'nodes={len(nodes)}  C3D10 elems={len(elems)}')
print(f'bore(FIX) nodes={len(fix_nodes)}   loaded faces={len(loaded)}')
print(f'flat face z={zface:.3f}   active band r in [{R_ACT_IN},{max(math.hypot(x,y) for _,x,y,z in nodes):.1f}]')
# through-thickness check
tvals=[z for _,_,_,z in nodes if abs(math.hypot(_,_)*0)==0]  # dummy
gmsh.finalize()
