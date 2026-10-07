"""
AXIAL MOTOR - multi-slice magnetic model, FULL SWEEP + convergence  BOX 3 of 4
=============================================================================
Reuses the working single slice and loops it over radius. FEMM is 2D and the
axial field is 3D, so we cut the active band into rings, unroll each into a flat
loaded linear model, solve it, and stack the per-slice axial force weighted by
the area each ring represents. Runs several slice counts to show convergence.

Per slice the geometry is the SAME shape that already worked at R=45mm; only the
tooth pitch and magnet pitch scale with that ring's radius (bigger ring = wider
spacing), and the slice depth = band width / number of slices.

Active radial band: 21 to 69 mm (from CAD).
Materials (exact names): N40, M-19 Steel, 18 AWG, Air.  Loaded: 432 A-turns.

RUN:  python axial_multislice.py     (FEMM installed; you drive it)
It prints total axial force for slice counts 5, 9, 15 -> watch it converge.
"""

from paths import *  # repo paths; also switches to the work folder
import os
import femm, math

R_IN, R_OUT = 21.0, 69.0          # active band (mm)
N_TEETH_TILE, N_POLES_TILE = 3, 2 # clean tile
TOOTH_H  = 14.33
MAGNET_T = 2.15
AIRGAP   = 1.80
AMP_TURNS = 432.0
SLICE_COUNTS = [5, 9, 15]

def build_and_solve_slice(R, depth):
    """Build the unrolled loaded slice at radius R (mm), return axial force (N)."""
    tooth_pitch  = 2*math.pi*R / 12.0     # 12 teeth around
    magnet_pitch = 2*math.pi*R / 8.0      # 8 poles around
    tooth_w  = tooth_pitch  * 0.55
    magnet_w = magnet_pitch * 0.80
    win_w = N_TEETH_TILE * tooth_pitch    # = 2 * magnet_pitch (clean tile)

    femm.newdocument(0)
    femm.mi_probdef(0,'millimeters','planar',1e-8,depth,30)
    for m in ('Air','N40','M-19 Steel','18 AWG'): femm.mi_getmaterial(m)

    def rect(x0,y0,x1,y1):
        femm.mi_drawline(x0,y0,x1,y0); femm.mi_drawline(x1,y0,x1,y1)
        femm.mi_drawline(x1,y1,x0,y1); femm.mi_drawline(x0,y1,x0,y0)
    def label(x,y,mat,circ='<None>',magdir=0,turns=0):
        femm.mi_addblocklabel(x,y); femm.mi_selectlabel(x,y)
        femm.mi_setblockprop(mat,1,0,circ,magdir,0,turns); femm.mi_clearselected()

    y0=0.0; y_bmag_t=MAGNET_T; y_tb=y_bmag_t+AIRGAP; y_tt=y_tb+TOOTH_H
    y_tmag_b=y_tt+AIRGAP; y_top=y_tmag_b+MAGNET_T
    rect(0,y0,win_w,y_top)

    for j in range(N_POLES_TILE):
        cx=(j+0.5)*magnet_pitch; mx0=cx-magnet_w/2; pol=90 if j%2==0 else 270
        rect(mx0,y0,mx0+magnet_w,y_bmag_t);      label(cx,y0+MAGNET_T/2,'N40','<None>',pol)
        rect(mx0,y_tmag_b,mx0+magnet_w,y_top);   label(cx,y_tmag_b+MAGNET_T/2,'N40','<None>',pol)

    tooth_cx=[(i+0.5)*tooth_pitch for i in range(N_TEETH_TILE)]
    for cx in tooth_cx:
        rect(cx-tooth_w/2,y_tb,cx+tooth_w/2,y_tt); label(cx,(y_tb+y_tt)/2,'M-19 Steel')

    femm.mi_addcircprop('A',AMP_TURNS,1); femm.mi_addcircprop('B',-AMP_TURNS,1)
    for (cl,cr,ph) in [(tooth_cx[0],tooth_cx[1],'A'),(tooth_cx[1],tooth_cx[2],'B')]:
        sx0=cl+tooth_w/2; sx1=cr-tooth_w/2
        rect(sx0,y_tb,sx1,y_tt); label((sx0+sx1)/2,(y_tb+y_tt)/2,'18 AWG',ph,0,1)

    label(min(1.0,win_w*0.02), y_bmag_t+AIRGAP*0.5, 'Air')

    femm.mi_addboundprop('A0',0,0,0,0,0,0,0,0,0)
    for (sx,sy) in [(win_w/2,y0),(win_w/2,y_top),(0,y_top/2),(win_w,y_top/2)]:
        femm.mi_selectsegment(sx,sy)
    femm.mi_setsegmentprop('A0',0,1,0,0); femm.mi_clearselected()

    femm.mi_saveas(os.path.join(WORK, '_axial_tmp.fem'))   # FEMM requires a saved file to solve
    femm.mi_analyze(1); femm.mi_loadsolution()
    for j in range(N_POLES_TILE):
        femm.mo_selectblock((j+0.5)*magnet_pitch, y_tmag_b+MAGNET_T/2)
    fy=femm.mo_blockintegral(19); femm.mo_clearblock()
    femm.mo_close()
    # scale the tile force up to the full ring: ring has 8 poles, tile has 2
    return fy * (8.0/N_POLES_TILE)

def main():
    femm.openfemm(0)
    print("radius-swept axial force (per-ring, scaled to full 8 poles):\n")
    for N in SLICE_COUNTS:
        edges=[R_IN+(R_OUT-R_IN)*k/N for k in range(N+1)]
        depth=(R_OUT-R_IN)/N
        total=0.0
        for k in range(N):
            Rmid=0.5*(edges[k]+edges[k+1])
            f=build_and_solve_slice(Rmid, depth)
            total+=f
        print(f"  {N:2d} slices -> total axial force = {total:8.2f} N")
    print("\nif the last two totals are close, the slice count has converged.")
    input("done. press Enter to close FEMM...")

if __name__=='__main__':
    main()
