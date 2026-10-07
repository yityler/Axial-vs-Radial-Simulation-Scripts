"""
Build the axial-disc harmonic (steady-state dynamics) deck + pressure wave.  BOX 4
Mirrors the radial pipeline: Step-1 *Frequency (storage), Step-2 *Steady state
dynamics + *Modal damping, load = per-face axial pressure wave on the disc face.

Axial force spectrum at 400 Hz (2x fe), de-aliased from axial_force_sweep.csv:
  disc order 0  : 16400 N/m^2   (uniform breathing - dominant drumhead driver)
  disc order 8  :  5347 N/m^2   (pole-count)
  disc order 12 :  5716 N/m^2   (slot-count)
Pressure p(theta) = [A0 + A8 cos(8 th) + A12 cos(12 th)]  (N/m^2 -> MPa via /1e6)
(phases unknown from magnitude-only FFT; cos with zero phase, first pass.)
"""
from paths import *  # repo paths; also switches to the work folder
import math, sys
SUF = sys.argv[1] if len(sys.argv)>1 else ''
A0, A8, A12 = 16400.0, 5347.0, 5716.0
PA_TO_MPA = 1e-6

faces = []
for ln in open(rf'disc_faces{SUF}.txt'):
    e,f,th,r = ln.split(',')
    faces.append((int(e),int(f),float(th)))

dl = []
for e,f,th in faces:
    p = (A0 + A8*math.cos(8*th) + A12*math.cos(12*th)) * PA_TO_MPA
    dl.append(f'{e}, P{f}, {p:.8e}')

deck = r'''*Heading
Axial rotor disc - harmonic (steady-state dynamics) at 400 Hz - Box 4
*Include, Input=disc_mesh{SUF}.inp
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
390, 410, 11, 1
*Modal damping
1, 20, 0.02
*Dload, Load case=1
'''
deck += '\n'.join(dl) + '\n'
deck += '''*Node file
U
*End step
'''
deck = deck.replace('{SUF}', SUF)
open(rf'disc_harm{SUF}.inp','w').write(deck)
print(f'wrote disc_harm{SUF}.inp : {len(dl)} loaded faces')
print(f'pressure range: min={min((A0+A8*math.cos(8*t)+A12*math.cos(12*t)) for _,_,t in faces):.0f}  max={max((A0+A8*math.cos(8*t)+A12*math.cos(12*t)) for _,_,t in faces):.0f} N/m^2')
