"""
RADIAL MOTOR - loaded magnetic run, force extraction  (BOX 1 of 4)
==================================================================
Fully wired version. Slot coordinates come from your CAD and were confirmed
against FEMM (rightmost slot at (38,0) matches CAD (37.5,0), no rotation offset).

WHAT IT DOES
  Energizes the 12-slot / 8-pole winding at the rated EV operating point,
  steps the rotor through one electrical period, and pulls the RADIAL air-gap
  Maxwell stress. It FFTs that force in space (around the gap) and time (rotor
  steps) to give the radial force harmonics by spatial order and frequency.
  Output: radial_loaded_force.csv  +  printed dominant harmonics.

SETTINGS  (locked by you, EV-appropriate; sources to add to Works Cited)
  current density   6 A/mm^2     -> continuous rating, liquid-cooled traction
                                    [CITE: motor-design text / EV motor paper]
  current angle     0 deg        -> max torque/amp for SURFACE-PM rotor
                                    [CITE: standard PMSM / MTPA reference]
  rated speed       3000 rpm     -> representative rated point (your choice,
                                    only used to convert orders -> Hz)
  amp-turns/slot    432          -> 72 mm^2 slot copper x 6 A/mm^2 (your CAD)
  winding (Emetor, 12s/8p, 2-layer): aB bC cA aB bC cA aB bC cA aB bC cA

HOW TO RUN
  pip install pyfemm numpy      (once)
  python radial_loaded_force_v2.py
  You drive FEMM; I cannot. If it errors, send me the traceback.
"""

from paths import *  # repo paths; also switches to the work folder
import femm, numpy as np, csv, os, math

# ---- paths ----------------------------------------------------------------
FEM_PATH = os.path.join(DATA, 'femm', 'radial_flux.FEM')
OUT_CSV  = os.path.join(DATA, 'radial_loaded_force.csv')

# ---- locked settings ------------------------------------------------------
AMP_TURNS      = 432.0     # A-turns per slot
CURRENT_ANGLE  = 0.0       # elec deg (surface PM -> pure q-axis)
RATED_RPM      = 3000.0
N_POLES        = 8
P              = N_POLES // 2

# ---- slot geometry (from CAD, confirmed vs FEMM) --------------------------
SLOT_R_MM   = 37.5
SLOT_ANGLES = [0,30,60,90,120,150,180,210,240,270,300,330]   # deg, slot 1..12

# ---- verified winding (Emetor, 12s/8p, SINGLE layer) ----------------------
# Layout string read off the Emetor diagram:  A a C c B b A a C c B b
# One phase per slot (matches your one-copper-region-per-slot FEMM geometry).
# lower = +, upper = -.  Verified: order-4 (8-pole) field dominates, each phase
# balanced (4 slots, net sign 0).  [CITE: Emetor winding calculator]
SLOT_LAYOUT = ["A","a","C","c","B","b","A","a","C","c","B","b"]

def slot_drive(ch):
    """Return (phase, sign) for this slot's single-layer coil."""
    return ch.upper(), (+1 if ch.islower() else -1)

SLOT_DRIVE = [slot_drive(c) for c in SLOT_LAYOUT]

# ---- sweep ----------------------------------------------------------------
N_STEPS        = 24                # finer temporal resolution -> trustworthy frequencies
MECH_PERIOD    = 360.0 / P          # one elec period in mech deg = 90
AIRGAP_R_MM    = 36.0               # just inside the slot radius; adjust if needed
N_CONTOUR_PTS  = 180               # points around the gap (plenty for spatial orders)

def phase_currents(elec_deg):
    th = math.radians(elec_deg + CURRENT_ANGLE)
    return (AMP_TURNS*math.cos(th),
            AMP_TURNS*math.cos(th-2*math.pi/3),
            AMP_TURNS*math.cos(th+2*math.pi/3))

def main():
    if not os.path.exists(FEM_PATH):
        raise FileNotFoundError(FEM_PATH)
    femm.openfemm(0)   # visible, so any solver dialog is clickable (not a hidden hang)
    print("FEMM opened, loading model...", flush=True)
    femm.opendocument(FEM_PATH)
    print("model loaded, wiring winding to slots...", flush=True)
    femm.mi_saveas(os.path.join(WORK, 'radial_flux_loaded.FEM'))

    # ensure 3 series circuits exist
    for name in ("A","B","C"):
        try: femm.mi_addcircprop(name, 0, 1)
        except Exception: pass

    # attach each slot's block label to its phase circuit with the right sign.
    # we select the block label nearest each slot centroid and set its circuit.
    for (ang, (phase, sign)) in zip(SLOT_ANGLES, SLOT_DRIVE):
        x = SLOT_R_MM*math.cos(math.radians(ang))
        y = SLOT_R_MM*math.sin(math.radians(ang))
        femm.mi_selectlabel(x, y)
        # turns carry the sign; magnitude 1 turn (amp-turns live in the current)
        femm.mi_setblockprop("18 AWG", 1, 0, phase, 0, 0, sign*1)
        femm.mi_clearselected()

    angles = np.linspace(0,360,N_CONTOUR_PTS,endpoint=False)
    Fmat = np.zeros((N_STEPS, N_CONTOUR_PTS))
    rows = []
    print("wiring done, starting rotor sweep...", flush=True)

    for k in range(N_STEPS):
        print(f"  step {k+1}/{N_STEPS}: setting currents...", flush=True)
        mech = MECH_PERIOD*k/N_STEPS
        elec = P*mech
        Ia,Ib,Ic = phase_currents(elec)
        femm.mi_setcurrent("A",Ia); femm.mi_setcurrent("B",Ib); femm.mi_setcurrent("C",Ic)

        femm.mi_selectgroup(1)                      # rotor = group 1 (as in cogging run)
        femm.mi_moverotate(0,0, MECH_PERIOD/N_STEPS if k>0 else 0)
        femm.mi_clearselected()

        print(f"  step {k+1}/{N_STEPS}: solving...", flush=True)
        femm.mi_analyze(1); femm.mi_loadsolution()   # 1 = solve without popup window
        print(f"  step {k+1}/{N_STEPS}: reading air-gap field...", flush=True)   # 1 = solve without popup window

        mu0 = 4*math.pi*1e-7
        for j,a in enumerate(angles):
            x = AIRGAP_R_MM*math.cos(math.radians(a))
            y = AIRGAP_R_MM*math.sin(math.radians(a))
            bx,by = femm.mo_getb(x,y)
            br =  bx*math.cos(math.radians(a)) + by*math.sin(math.radians(a))
            bt = -bx*math.sin(math.radians(a)) + by*math.cos(math.radians(a))
            sr = (br*br - bt*bt)/(2*mu0)            # radial Maxwell stress N/m^2
            Fmat[k,j] = sr
            rows.append([k,elec,a,br,bt,sr])

    femm.closefemm()

    with open(OUT_CSV,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["step","elec_deg","angle_deg","Br","Bt","sigma_r_Npm2"])
        w.writerows(rows)
    print("wrote", OUT_CSV)

    F = np.fft.rfft2(Fmat); mag = np.abs(F)/Fmat.size
    fe = RATED_RPM/60.0 * P
    flat = [(mag[n,r],r,n) for n in range(mag.shape[0]) for r in range(mag.shape[1])]
    flat.sort(reverse=True)
    print("\nDominant radial-force harmonics (spatial order r, temporal order n):")
    for amp,r,n in flat[1:11]:
        print(f"  order r={r:3d}  temporal n={n:2d}  ~{n*fe:8.1f} Hz   amp {amp:.3e} N/m^2")

if __name__ == "__main__":
    main()
