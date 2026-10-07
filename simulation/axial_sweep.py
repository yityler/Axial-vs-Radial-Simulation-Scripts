"""
AXIAL MOTOR - rotor sweep, axial force pattern + harmonics   BOX 3 -> load for BOX 4
===================================================================================
Extends the working slice to produce the axial force DISTRIBUTION (not just the
total), so the structural run has a real force wave to apply. Uses the middle
slice (R=45 mm). At each rotor position it samples the AXIAL Maxwell stress along
the top air gap vs position, then FFTs in space (around the disc) and time
(rotor position) to give the axial force harmonics.

Output: axial_force_sweep.csv  +  printed dominant axial harmonics.

FIRST-PASS SIMPLIFICATIONS (flagged honestly):
  - single (middle) slice used for the force SHAPE; total magnitude already known
    (~242 N converged from the multi-slice run)
  - phase currents held fixed during the sweep (rotor-only). Advancing current
    synchronously is a refinement.
  - magnet arrangement + estimated widths + A=0 boundary carried from the slice,
    still to be verified against a reference / Dr. Anwar.

RUN:  python axial_sweep.py   (run twice if the first errors on materials)
"""
from paths import *  # repo paths; also switches to the work folder
import os
import femm, math

R_MID=45.0
TOOTH_PITCH=2*math.pi*R_MID/12.0
MAGNET_PITCH=2*math.pi*R_MID/8.0
TOOTH_W=TOOTH_PITCH*0.55
MAGNET_W=MAGNET_PITCH*0.80
TOOTH_H=14.33; MAGNET_T=2.15; AIRGAP=1.80; DEPTH=16.0
AMP_TURNS=432.0
N_TEETH=6                      # 2 tiles wide, central tile is the clean sample zone
W=N_TEETH*TOOTH_PITCH
N_STEPS=24                     # rotor positions over one electrical period
SWEEP=2*MAGNET_PITCH           # one electrical period = 2 pole pitches
N_SAMP=120                     # points sampled along the gap (central window)
MU0=4*math.pi*1e-7

def rect(x0,y0,x1,y1):
    femm.mi_drawline(x0,y0,x1,y0); femm.mi_drawline(x1,y0,x1,y1)
    femm.mi_drawline(x1,y1,x0,y1); femm.mi_drawline(x0,y1,x0,y0)
def label(x,y,mat,circ='<None>',magdir=0,turns=0):
    femm.mi_addblocklabel(x,y); femm.mi_selectlabel(x,y)
    femm.mi_setblockprop(mat,1,0,circ,magdir,0,turns); femm.mi_clearselected()

def build(offset):
    femm.newdocument(0)
    femm.mi_probdef(0,'millimeters','planar',1e-8,DEPTH,30)
    for m in ('Air','N40','M-19 Steel','18 AWG'): femm.mi_getmaterial(m)
    y0=0.0; y_bmt=MAGNET_T; y_tb=y_bmt+AIRGAP; y_tt=y_tb+TOOTH_H
    y_tmb=y_tt+AIRGAP; y_top=y_tmb+MAGNET_T
    rect(0,y0,W,y_top)
    # magnets tiled across the box, shifted by rotor offset, polarity by pole parity
    m=-1
    while (m+0.5)*MAGNET_PITCH+offset < W+MAGNET_PITCH:
        cx=(m+0.5)*MAGNET_PITCH+offset
        if MAGNET_W/2 < cx < W-MAGNET_W/2:
            pol=90 if m%2==0 else 270
            rect(cx-MAGNET_W/2,y0,cx+MAGNET_W/2,y_bmt);   label(cx,y0+MAGNET_T/2,'N40','<None>',pol)
            rect(cx-MAGNET_W/2,y_tmb,cx+MAGNET_W/2,y_top);label(cx,y_tmb+MAGNET_T/2,'N40','<None>',pol)
        m+=1
    # stator teeth (fixed)
    tcx=[(i+0.5)*TOOTH_PITCH for i in range(N_TEETH)]
    for cx in tcx:
        rect(cx-TOOTH_W/2,y_tb,cx+TOOTH_W/2,y_tt); label(cx,(y_tb+y_tt)/2,'M-19 Steel')
    # coils in interior slots, fixed currents (loaded)
    femm.mi_addcircprop('A',AMP_TURNS,1); femm.mi_addcircprop('B',-AMP_TURNS,1)
    for i in range(N_TEETH-1):
        sx0=tcx[i]+TOOTH_W/2; sx1=tcx[i+1]-TOOTH_W/2
        ph='A' if i%2==0 else 'B'
        rect(sx0,y_tb,sx1,y_tt); label((sx0+sx1)/2,(y_tb+y_tt)/2,'18 AWG',ph,0,1)
    label(min(1.0,W*0.01), y_bmt+AIRGAP*0.5,'Air')
    femm.mi_addboundprop('A0',0,0,0,0,0,0,0,0,0)
    for (sx,sy) in [(W/2,y0),(W/2,y_top),(0,y_top/2),(W,y_top/2)]: femm.mi_selectsegment(sx,sy)
    femm.mi_setsegmentprop('A0',0,1,0,0); femm.mi_clearselected()
    femm.mi_saveas(os.path.join(WORK, '_axial_sweep_tmp.fem'))
    femm.mi_analyze(1); femm.mi_loadsolution()
    return y_tt+AIRGAP*0.5   # y of top-gap sampling line

def main():
    import numpy as np, csv
    femm.openfemm(0)
    x_lo=W*0.25; x_hi=W*0.75           # central window (edge effects excluded)
    xs=[x_lo+(x_hi-x_lo)*k/N_SAMP for k in range(N_SAMP)]
    Fmat=np.zeros((N_STEPS,N_SAMP)); rows=[]
    for s in range(N_STEPS):
        off=SWEEP*s/N_STEPS
        ygap=build(off)
        for j,x in enumerate(xs):
            bx,by=femm.mo_getb(x,ygap)
            sig=(by*by-bx*bx)/(2*MU0)   # axial (y) normal Maxwell stress
            Fmat[s,j]=sig; rows.append([s,off,x,bx,by,sig])
        femm.mo_close()
    with open(os.path.join(DATA, 'axial_force_sweep.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['step','offset_mm','x_mm','Bx','By','sigma_axial_Npm2']); w.writerows(rows)
    print('wrote axial_force_sweep.csv')
    F=np.fft.rfft2(Fmat); mag=np.abs(F)/Fmat.size
    fe=200.0
    flat=[(mag[n,r],r,n) for n in range(mag.shape[0]) for r in range(mag.shape[1])]
    flat.sort(reverse=True)
    print('\nDominant AXIAL force harmonics (window order r, temporal n, freq):')
    print('(disc spatial order = window order x 4, since window = 1 tile = 2 of 8 poles)')
    c=0
    for amp,r,n in flat:
        if n==0: continue
        if n>N_STEPS//2: continue   # axis-0 is a full FFT; bins n and N_STEPS-n are conjugate mirrors
        print(f'  window r={r:2d} (disc order {r*4:2d}) | n={n:2d} | {n*fe:6.0f} Hz | amp {amp:.3e} N/m^2')
        c+=1
        if c>=8: break
    input('done. press Enter to close FEMM...')

if __name__=='__main__':
    main()
