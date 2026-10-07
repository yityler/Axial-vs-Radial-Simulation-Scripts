"""
AXIAL MOTOR - multi-slice magnetic model, SLICE 1 (middle radius)   BOX 3 of 4
=============================================================================
Rebuilt around the CORRECT repeating unit. For 12 slots / 8 poles the smallest
clean tile is 3 teeth + 2 poles (~70.7 mm), so whole teeth and whole magnets
fit with nothing cut off at the window edges.

Unrolled ring at R = 45 mm: two rotor magnet rows (top+bottom), stator teeth in
the middle, air gaps between. Solves and reports the AXIAL force (Fy) on the
teeth, which drives axial (drumhead) vibration. This is slice 1 of 3; once it
looks right we loop it over R = 27, 45, 63 mm.

REAL DIMENSIONS at R = 45 mm (from your Fusion CAD):
  tooth pitch 23.56 mm (x3)   tooth height 14.33 mm
  magnet pitch 35.34 mm (x2)  magnet thick 2.15 mm   air gap ~1.8 mm/side
  materials (exact model names): N40, M-19 Steel, 18 AWG, Air

FLAGGED ASSUMPTIONS (confirm later / with Dr. Anwar):
  - tooth width = 55% of pitch, magnet width = 80% of pole pitch (typical arcs)
  - magnet rows arranged N-S facing across the gap (axial flux through teeth)
  - outer box held at A=0 (approximates rotor back-iron surface) for this pass
  - slice depth (out-of-plane) = 16 mm (48 mm active band / 3 slices)

RUN:  python axial_slice1.py    (FEMM must be installed; you drive it)
"""

from paths import *  # repo paths; also switches to the work folder
import os
import femm, math

R_MID        = 45.0
TOOTH_PITCH  = 23.56
MAGNET_PITCH = 35.34
TOOTH_H      = 14.33
MAGNET_T     = 2.15
AIRGAP       = 1.80
TOOTH_W      = TOOTH_PITCH * 0.55
MAGNET_W     = MAGNET_PITCH * 0.80
SLICE_DEPTH  = 16.0

N_TEETH   = 3                     # clean tile: 3 teeth
N_POLES   = 2                     # clean tile: 2 poles
WIN_W     = N_TEETH * TOOTH_PITCH # = 2 * MAGNET_PITCH = 70.68 mm

def rect(x0,y0,x1,y1):
    femm.mi_drawline(x0,y0,x1,y0); femm.mi_drawline(x1,y0,x1,y1)
    femm.mi_drawline(x1,y1,x0,y1); femm.mi_drawline(x0,y1,x0,y0)

def label(x,y,mat,magdir=0):
    femm.mi_addblocklabel(x,y); femm.mi_selectlabel(x,y)
    femm.mi_setblockprop(mat,1,0,'<None>',magdir,0,0); femm.mi_clearselected()

def main():
    femm.openfemm(0)
    femm.newdocument(0)
    femm.mi_probdef(0,'millimeters','planar',1e-8,SLICE_DEPTH,30)

    for m in ('Air','N40','M-19 Steel','18 AWG'):
        femm.mi_getmaterial(m)

    # vertical layout: bottom magnets / gap / teeth / gap / top magnets
    y0        = 0.0
    y_bmag_t  = y0 + MAGNET_T
    y_teeth_b = y_bmag_t + AIRGAP
    y_teeth_t = y_teeth_b + TOOTH_H
    y_tmag_b  = y_teeth_t + AIRGAP
    y_top     = y_tmag_b + MAGNET_T

    rect(0,y0,WIN_W,y_top)                     # outer box

    # --- magnet rows (2 poles), alternating polarity, N-S facing across gap ---
    for j in range(N_POLES):
        cx = (j+0.5)*MAGNET_PITCH
        mx0 = cx - MAGNET_W/2
        pol = 90 if j%2==0 else 270            # alternate poles
        # bottom row
        rect(mx0,y0,mx0+MAGNET_W,y_bmag_t)
        label(cx, y0+MAGNET_T/2, 'N40', pol)
        # top row (same angle -> presents opposite face to the gap = N-S facing)
        rect(mx0,y_tmag_b,mx0+MAGNET_W,y_top)
        label(cx, y_tmag_b+MAGNET_T/2, 'N40', pol)

    # --- teeth (3) ---
    tooth_x = []
    for i in range(N_TEETH):
        cx = (i+0.5)*TOOTH_PITCH
        tooth_x.append(cx)
        tx0 = cx - TOOTH_W/2
        rect(tx0,y_teeth_b,tx0+TOOTH_W,y_teeth_t)
        label(cx, (y_teeth_b+y_teeth_t)/2, 'M-19 Steel', 0)

    # --- coils drawn as enclosed rectangles in the two interior slots ---
    # a slot is the gap between two adjacent teeth; drawing a bounded rectangle
    # there makes it its own region so copper does not leak into the air.
    femm.mi_addcircprop('A',  432.0, 1)
    femm.mi_addcircprop('B', -432.0, 1)
    interior = [(tooth_x[0], tooth_x[1], 'A'),
                (tooth_x[1], tooth_x[2], 'B')]
    for (cl, cr, ph) in interior:
        sx0 = cl + TOOTH_W/2          # right edge of left tooth
        sx1 = cr - TOOTH_W/2          # left edge of right tooth
        rect(sx0, y_teeth_b, sx1, y_teeth_t)
        cxm = (sx0+sx1)/2
        femm.mi_addblocklabel(cxm, (y_teeth_b+y_teeth_t)/2)
        femm.mi_selectlabel(cxm, (y_teeth_b+y_teeth_t)/2)
        femm.mi_setblockprop('18 AWG', 1, 0, ph, 0, 0, 1)
        femm.mi_clearselected()

    # --- air: one label in the lower gap (all gap/slot air is connected) ---
    label(2.0, y_bmag_t+AIRGAP*0.5, 'Air', 0)

    # --- outer boundary A=0 (bounds the problem so it solves) ---
    femm.mi_addboundprop('A0',0,0,0,0,0,0,0,0,0)   # prop type 0 = prescribed A=0
    for (sx,sy) in [(WIN_W/2,y0),(WIN_W/2,y_top),(0,y_top/2),(WIN_W,y_top/2)]:
        femm.mi_selectsegment(sx,sy)
    femm.mi_setsegmentprop('A0',0,1,0,0)
    femm.mi_clearselected()

    femm.mi_zoomnatural()
    femm.mi_saveas(os.path.join(WORK, 'axial_slice1.fem'))
    femm.mi_analyze(1)
    femm.mi_loadsolution()

    # axial force per side: integrate on the TOP magnet row (force it feels from
    # the teeth across the top gap). Net force on teeth cancels between the two
    # rotors by design, so we measure one interface instead.
    for j in range(N_POLES):
        cx = (j+0.5)*MAGNET_PITCH
        femm.mo_selectblock(cx, y_tmag_b + MAGNET_T/2)
    fy = femm.mo_blockintegral(19)   # 19 = y force (weighted stress tensor)
    femm.mo_clearblock()
    print('SLICE 1 (R=45mm) axial force across TOP gap  Fy = %.4f N' % fy)
    print('check FEMM: 2 magnet poles top+bottom, 3 teeth + 2 coils middle.')
    input('FEMM left open. Look, then press Enter here to close...')

if __name__ == '__main__':
    main()
