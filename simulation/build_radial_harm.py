"""
Build the radial-stator harmonic deck (gmsh+ccx pipeline).  Same-tooling rebuild.
Radial pressure wave on the bore: p(theta) = A4 cos(4 th) + A16 cos(16 th)
Amplitudes from radial_loaded_force.csv FFT (n=1, 200 Hz):
  order 4  = 3452 N/m^2 ,  order 16 = 4117 N/m^2
Outer rim fixed. Sweep 190-210 Hz, 2% modal damping, 20 modes.
Usage: python build_radial_harm.py [suffix]
"""
from paths import *  # repo paths; also switches to the work folder
import math, sys
SUF = sys.argv[1] if len(sys.argv)>1 else ''
A4, A16 = 3452.0, 4117.0
PA_TO_MPA = 1e-6

faces=[]
for ln in open(rf'radial_faces{SUF}.txt'):
    e,f,th,r = ln.split(','); faces.append((int(e),int(f),float(th)))

dl=[]
for e,f,th in faces:
    p = (A4*math.cos(4*th) + A16*math.cos(16*th))*PA_TO_MPA
    dl.append(f'{e}, P{f}, {p:.8e}')

deck = r'''*Heading
Radial stator - harmonic (steady-state dynamics) at 200 Hz - gmsh+ccx rebuild
*Include, Input=radial_mesh{SUF}.inp
*Material, Name=Steel
*Density
7.85E-09
*Elastic
200000, 0.3
*Solid section, Elset=SOLID, Material=Steel
*Boundary
FIX, 1, 3
*Step
*Frequency, Solver=Pardiso, Storage=Yes
20
*Node file
U
*End step
*Step
*Steady state dynamics, Solver=Pardiso
190, 210, 11, 1
*Modal damping
1, 20, 0.02
*Dload, Load case=1
'''
deck = deck.replace('{SUF}', SUF)
deck += '\n'.join(dl) + '\n'
deck += '''*Node file
U
*End step
'''
open(rf'radial_harm{SUF}.inp','w').write(deck)
print(f'wrote radial_harm{SUF}.inp : {len(dl)} loaded faces')
