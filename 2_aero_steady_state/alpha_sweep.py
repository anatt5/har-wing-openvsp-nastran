"""
Steady-state alpha sweep: aerodynamic derivatives extraction (VSPAERO).

Builds the aircraft geometry (main wing, horizontal tail, vertical tail,
control surfaces) at the fixed mesh density (d=50 panels/m) selected by the
separate VSPAERO mesh-convergence study, then sweeps angle of attack from
-4 deg to 12 deg (linear VLM, no stall model, undeflected controls) at a
fixed centre of gravity. From the resulting CL/CD/CMy vs alpha data it
fits:
  - the lift curve   (CL = CL_alpha*alpha + CL0)
  - the pitching-moment curve (Cm = Cm_alpha*alpha + Cm0)
  - the parabolic drag polar (CD = CD0 + k*CL^2)
and derives CL_alpha, alpha_L0, CL0, CD0, k, (L/D)max, Cm_alpha, the static
margin and the neutral point about the reference CG.

Outputs:
  AC.vsp3                 -- clean geometry (fixed d=50 mesh)
  log_vsp_trim.txt         -- full console log
  lift_curve.png, drag_polar.png, pitching_moment.png, lift_to_drag.png

Run:  python alpha_sweep.py
"""
import os
import sys
import openvsp as vsp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
AIRFOIL_DIR = os.path.join(REPO_ROOT, "airfoils")

sys.path.insert(0, REPO_ROOT)
from utils.airfoil_utils import assign_airfoil_to_component

vsp.VSPRenew()
vsp.VSPCheckSetup()


########################################################################
# 2) MAIN WING
########################################################################
wing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(wing_id, "main wing")
vsp.SetParmVal(wing_id, "Tess_W", "Shape", 49)

span_mw = 1.6
AR_mw = 6.4
taper_mw = 0.66667

vsp.SetDriverGroup(wing_id, 1,
                    vsp.SPAN_WSECT_DRIVER,
                    vsp.AR_WSECT_DRIVER,
                    vsp.TAPER_WSECT_DRIVER)
vsp.Update()

sweep_angle = 28
vsp.SetParmVal(vsp.GetParm(wing_id, "Span", "XSec_1"), span_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Aspect", "XSec_1"), AR_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Taper", "XSec_1"), taper_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Sweep", "XSec_1"), sweep_angle)
vsp.Update()

# wing position
vsp.SetParmVal(vsp.GetParm(wing_id, "X_Rel_Location", "XForm"), 0.7)
vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Location", "XForm"), 0)
vsp.SetParmVal(vsp.GetParm(wing_id, "Z_Rel_Location", "XForm"), 0)
vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Rotation", "XForm"), 0)
vsp.Update()

# Airfoil: Clark Y
assign_airfoil_to_component(wing_id, "clark Y.dat", AIRFOIL_DIR)
xsec_surf_id_mwing = vsp.GetXSecSurf(wing_id, 0)

xsec_id_1_mwing = vsp.GetXSec(xsec_surf_id_mwing, 0)
parm_id_1_wing = vsp.GetXSecParm(xsec_id_1_mwing, "SectTess_U")
vsp.SetParmVal(parm_id_1_wing, 86.0)
dihedral_id_1 = vsp.GetXSecParm(xsec_id_1_mwing, "Dihedral")
vsp.SetParmVal(dihedral_id_1, 6.0)

xsec_id_2_mwing = vsp.GetXSec(xsec_surf_id_mwing, 1)
parm_id_2_mwing = vsp.GetXSecParm(xsec_id_2_mwing, "SectTess_U")
vsp.SetParmVal(parm_id_2_mwing, 86.0)
dihedral_id_2 = vsp.GetXSecParm(xsec_id_2_mwing, "Dihedral")
vsp.SetParmVal(dihedral_id_2, 6.0)

vsp.SetParmVal(wing_id, "Density", "Mass_Props", 40)
vsp.Update()

########################################################################
# 3) HORIZONTAL TAIL (HTP)
########################################################################
htail_id = vsp.AddGeom("WING")
vsp.SetGeomName(htail_id, "HTP")
vsp.SetParmVal(htail_id, "Tess_W", "Shape", 53)

span_htp = 0.38
rootc_htp = 0.35
taper_htp = 0.57143

vsp.SetDriverGroup(htail_id, 1,
                    vsp.SPAN_WSECT_DRIVER,
                    vsp.TAPER_WSECT_DRIVER,
                    vsp.ROOTC_WSECT_DRIVER)
vsp.Update()

vsp.SetParmVal(vsp.GetParm(htail_id, "Span", "XSec_1"), span_htp)
vsp.SetParmVal(vsp.GetParm(htail_id, "Taper", "XSec_1"), taper_htp)
vsp.SetParmVal(vsp.GetParm(htail_id, "Root_Chord", "XSec_1"), rootc_htp)
vsp.SetParmVal(vsp.GetParm(htail_id, "Sweep", "XSec_1"), sweep_angle)

# HTP position (above the VTP)
vsp.SetParmVal(vsp.GetParm(htail_id, "X_Rel_Location", "XForm"), 2.55)
vsp.SetParmVal(vsp.GetParm(htail_id, "Z_Rel_Location", "XForm"), 0.615)
vsp.SetParmVal(vsp.GetParm(htail_id, "Y_Rel_Location", "XForm"), 0.0)
vsp.SetParmVal(vsp.GetParm(htail_id, "Y_Rel_Rotation", "XForm"), -3)
vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Location", "XForm"), 0)
vsp.Update()

# Airfoil: NACA 0012
assign_airfoil_to_component(htail_id, "naca 0012.dat", AIRFOIL_DIR)
xsec_surf_id_htp = vsp.GetXSecSurf(htail_id, 0)

xsec_id_1_htp = vsp.GetXSec(xsec_surf_id_htp, 0)
parm_id_1_htp = vsp.GetXSecParm(xsec_id_1_htp, "SectTess_U")
vsp.SetParmVal(parm_id_1_htp, 24.0)

xsec_id_2_htp = vsp.GetXSec(xsec_surf_id_htp, 1)
parm_id_2_htp = vsp.GetXSecParm(xsec_id_2_htp, "SectTess_U")
vsp.SetParmVal(parm_id_2_htp, 24.0)
vsp.SetParmVal(htail_id, "Density", "Mass_Props", 40)
xsec_surf_id = vsp.GetXSecSurf(htail_id, 0)
xsec_1 = vsp.GetXSec(xsec_surf_id, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1, "ThickChord"), 0.08)

########################################################################
# 4) VERTICAL TAIL (VTP) -- two sections (dorsal fin + main fin)
########################################################################
vtp_id = vsp.AddGeom("WING")
vsp.SetGeomName(vtp_id, "VTP")
vsp.SetParmVal(vtp_id, "Tess_W", "Shape", 80)

vsp.InsertXSec(vtp_id, 1, vsp.XS_FOUR_SERIES)
vsp.Update()

# Section 1 (dorsal fin, root)
vsp.SetDriverGroup(vtp_id, 1, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
vsp.Update()

span_vtp = 0.2
rootc_vtp = 0.5
taper_vtp = 0.68966
sweep_angle1 = 50

vsp.SetParmVal(vtp_id, "Span", "XSec_1", span_vtp)
vsp.SetParmVal(vtp_id, "Taper", "XSec_1", taper_vtp)
vsp.SetParmVal(vtp_id, "Root_Chord", "XSec_1", rootc_vtp)
vsp.SetParmVal(vtp_id, "Sweep", "XSec_1", sweep_angle1)

# Section 2 (main fin, tip)
vsp.SetDriverGroup(vtp_id, 2, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
vsp.Update()  # mandatory before changing section 2's drivers

span_vtp2 = 0.23
taper_vtp2 = 0.81034
rootc_vtp2 = 0.34483
sweep_angle2 = 35

vsp.SetParmVal(vtp_id, "Span", "XSec_2", span_vtp2)
vsp.SetParmVal(vtp_id, "Taper", "XSec_2", taper_vtp2)
vsp.SetParmVal(vtp_id, "Root_Chord", "XSec_2", rootc_vtp2)
vsp.SetParmVal(vtp_id, "Sweep", "XSec_2", sweep_angle2)

vsp.SetParmVal(vtp_id, "X_Rel_Rotation", "XForm", 90)
vsp.SetParmVal(vtp_id, "Sym_Planar_Flag", "Sym", 0)
vsp.SetParmVal(vtp_id, "X_Rel_Location", "XForm", 2.15)
vsp.SetParmVal(vtp_id, "Z_Rel_Location", "XForm", 0.17)
vsp.Update()

# Airfoil: NACA 0008
assign_airfoil_to_component(vtp_id, "naca 0008.dat", AIRFOIL_DIR)
xsec_surf_id = vsp.GetXSecSurf(vtp_id, 0)

# root (dorsal fin base)
xsec_0 = vsp.GetXSec(xsec_surf_id, 0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_0, "ThickChord"), 0.08)
vsp.SetParmVal(vsp.GetXSecParm(xsec_0, "SectTess_U"), 26.0)

# junction between dorsal fin and main fin
xsec_1 = vsp.GetXSec(xsec_surf_id, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1, "ThickChord"), 0.12)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1, "SectTess_U"), 26.0)

# tip
xsec_2 = vsp.GetXSec(xsec_surf_id, 2)
vsp.SetParmVal(vsp.GetXSecParm(xsec_2, "ThickChord"), 0.08)
vsp.SetParmVal(vsp.GetXSecParm(xsec_2, "SectTess_U"), 26.0)

vsp.Update()
vsp.SetParmVal(vtp_id, "Density", "Mass_Props", 40)

########################################################################
# 5) CLOSE THE GEOMETRY (round end caps)
########################################################################
vsp.SetParmVal(wing_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(wing_id, "CapUMaxOption", "EndCap", 1.0)
vsp.SetParmVal(htail_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(htail_id, "CapUMaxOption", "EndCap", 1.0)
vsp.SetParmVal(vtp_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(vtp_id, "CapUMaxOption", "EndCap", 1.0)

########################################################################
# 6) VSPAERO CONTROL SURFACES
########################################################################
aileron_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
vsp.SetParmVal(vsp.GetParm(aileron_id, "UStart", "SS_Control"), 0.55)
vsp.SetParmVal(vsp.GetParm(aileron_id, "UEnd", "SS_Control"), 0.65)
vsp.SetParmVal(vsp.GetParm(aileron_id, "Length_C_Start", "SS_Control"), 0.2)

flap_id = vsp.AddSubSurf(wing_id, vsp.SS_CONTROL)
vsp.SetParmVal(vsp.GetParm(flap_id, "UStart", "SS_Control"), 0.38)
vsp.SetParmVal(vsp.GetParm(flap_id, "UEnd", "SS_Control"), 0.55)
vsp.SetParmVal(vsp.GetParm(flap_id, "Length_C_Start", "SS_Control"), 0.2)

rudder_id = vsp.AddSubSurf(vtp_id, vsp.SS_CONTROL, 0)
vsp.SetParmVal(vsp.GetParm(rudder_id, "UStart", "SS_Control"), 0.3)
vsp.SetParmVal(vsp.GetParm(rudder_id, "UEnd", "SS_Control"), 0.7)
vsp.SetParmVal(vsp.GetParm(rudder_id, "Length_C_Start", "SS_Control"), 0.3)

elevator_id = vsp.AddSubSurf(htail_id, vsp.SS_CONTROL, 0)
vsp.SetParmVal(vsp.GetParm(elevator_id, "UStart", "SS_Control"), 0.35)
vsp.SetParmVal(vsp.GetParm(elevator_id, "UEnd", "SS_Control"), 0.63)
vsp.SetParmVal(vsp.GetParm(elevator_id, "Length_C_Start", "SS_Control"), 0.3)


def set_surface_gain_by_index(group_index, target_val_0, target_val_1):
    """Set the (symmetric/antisymmetric) gain pair for a control-surface
    group. Every group has exactly 2 gain parameters, in creation order."""
    all_settings = vsp.FindContainer("VSPAEROSettings", 0)
    all_parms = vsp.FindContainerParmIDs(all_settings)

    gain_pids = [pid for pid in all_parms if "Gain" in vsp.GetParmName(pid)]

    idx0 = group_index * 2
    idx1 = group_index * 2 + 1
    if idx1 >= len(gain_pids):
        print(f"ERROR: group {group_index} does not exist")
        return

    vsp.SetParmVal(gain_pids[idx0], target_val_0)
    vsp.SetParmVal(gain_pids[idx1], target_val_1)
    vsp.Update()


# ailerons: differential (roll)
g_ailerons_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Ailerons_Group", g_ailerons_index)
vsp.AddSelectedToCSGroup([1, 2], g_ailerons_index)
set_surface_gain_by_index(0, 1, 1)

# flaps: symmetric
g_flaps_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Flaps_Group", g_flaps_index)
vsp.AddSelectedToCSGroup([3, 4], g_flaps_index)
set_surface_gain_by_index(1, 1.0, -1.0)

# elevators: symmetric, undeflected for this clean-configuration sweep
g_elevators_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Elevators_Group", g_elevators_index)
vsp.AddSelectedToCSGroup([5, 6], g_elevators_index)
set_surface_gain_by_index(2, 1.0, -1.0)

container_id = vsp.FindContainer("VSPAEROSettings", 0)
group_name = f"ControlSurfaceGroup_{g_elevators_index}"
parm_id = vsp.FindParm(container_id, "DeflectionAngle", group_name)
elevator_deflection = 0.0   # clean configuration: control surfaces undeflected
vsp.SetParmVal(parm_id, elevator_deflection)
vsp.Update()

# rudder
g_rudder_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Rudder_Group", g_rudder_index)
vsp.AddSelectedToCSGroup([7], g_rudder_index)
vsp.Update()

########################################################################
# 7) FIXED MESH: d = 50 panels/m (selected by the mesh-convergence study)
########################################################################
# Uniform mesh chosen as grid-independent:
#   main wing : Num W = 25, Num U = 81
#   HTP       : Num W = 29, Num U = 20
#   VTP       : Num W = 37, Num U = 12 (both sections)
# Num W = Tess_W ("Shape"); Num U = SectTess_U (set on every xsec to be safe).
for geom_id, num_w, num_u in ((wing_id, 25, 81), (htail_id, 29, 20), (vtp_id, 37, 12)):
    vsp.SetParmVal(geom_id, "Tess_W", "Shape", num_w)
    _surf = vsp.GetXSecSurf(geom_id, 0)
    for _i in range(vsp.GetNumXSec(_surf)):
        vsp.SetParmVal(vsp.GetXSecParm(vsp.GetXSec(_surf, _i), "SectTess_U"), float(num_u))
vsp.Update()

# save clean geometry (no analysis artifacts)
vsp.WriteVSPFile("AC.vsp3")
print("plane geometry created: AC.vsp3")

########################################################################
# 8) STEADY-STATE ALPHA SWEEP
########################################################################
import sys


class Logger(object):
    """Duplicate every print() to a log file as well as the console."""

    def __init__(self, filename="log_vsp_trim.txt"):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()


sys.stdout = Logger("log_vsp_trim.txt")

# ======================== reference constants ========================
total_mass_kg = 20.0        # total mass (kg) -- TU-Flex (~20 kg, papers)
freestream_velocity = 30.0  # Vinf (m/s) -- TU-Flex cruise speed Vc (papers)
wing_area = 0.8              # Sref (m2)
wingspan = 3.2                # bref (m)

# centre of gravity
X_cg = 1.23   # reference CG (near the neutral point)
Y_cg = 0.00
Z_cg = 0.023

# air properties (ISA, sea level)
air_density = 1.225
speed_of_sound = 340.29
dynamic_viscosity = 1.789e-5
gravity = 9.80665

# elevator deflection for the sweep (0 = clean configuration)
elevator_deflection = 0.0

mean_chord = 0.25321   # MAC from OpenVSP (consistent with the mesh convergence cref)
mach_start = freestream_velocity / speed_of_sound
re_cref = (air_density * freestream_velocity * mean_chord) / dynamic_viscosity
g_elevators_index = 2

print("=" * 60)
print("ALPHA SWEEP INPUT CONFIGURATION:")
print(f"-> Freestream velocity:  {freestream_velocity:.1f} m/s")
print(f"-> CG position (X_cg):   {X_cg:.3f} m")
print(f"-> Elevator deflection:  {elevator_deflection:.4f} deg")
print(f"-> Mach number:          {mach_start:.4f}")
print(f"-> Reynolds (ReCref):    {re_cref:.1f}")
print("=" * 60)


def get_cl_cd_cm(filename="Aircraft_Stability.csv"):
    """Read the converged CLtot, CDtot and CMytot (last column) from the results CSV."""
    cl = cd = cm = None
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith('CLtot,'):
                    cl = float(s.split(',')[-1])
                elif s.startswith('CDtot,'):
                    cd = float(s.split(',')[-1])
                elif s.startswith('CMytot,'):
                    cm = float(s.split(',')[-1])
    return cl, cd, cm


# ======================== VSPAERO inputs ========================
analysis_name = "VSPAEROSweep"
vsp.SetAnalysisInputDefaults(analysis_name)

vsp.SetDoubleAnalysisInput(analysis_name, "Sref", [wing_area], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "bref", [wingspan], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "cref", [mean_chord], 0)

vsp.SetDoubleAnalysisInput(analysis_name, "Xcg", [X_cg], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Ycg", [Y_cg], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Zcg", [Z_cg], 0)

vsp.SetDoubleAnalysisInput(analysis_name, "Rho", [air_density], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", [freestream_velocity], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "MachStart", [mach_start], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "ReCref", [re_cref], 0)
vsp.SetIntAnalysisInput(analysis_name, "Symmetry", [0], 0)

# wake and solver settings (free relaxed wake)
vsp.SetIntAnalysisInput(analysis_name, "WakeNumIter", [5], 0)
vsp.SetIntAnalysisInput(analysis_name, "NumWakeNodes", [64], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "WakeRelax", [0.5], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "CoreSizeFactor", [1.0], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "ThinGeomSet", [0], 0)   # all surfaces thin -> VLM (matches mesh convergence)
vsp.SetDoubleAnalysisInput(analysis_name, "GeomSet", [0], 0)

vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [1], 0)
vsp.SetIntAnalysisInput(analysis_name, "BetaNpts", [1], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "BetaStart", [0.0], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "BetaEnd", [0.0], 0)

# locate the elevator parameter and apply its deflection
container_id = vsp.FindContainer("VSPAEROSettings", 0)
group_name = f"ControlSurfaceGroup_{g_elevators_index}"
parm_id = vsp.FindParm(container_id, "DeflectionAngle", group_name)
vsp.SetParmVal(parm_id, elevator_deflection)
vsp.Update()

# ======================== alpha sweep ========================
print("\n" + "=" * 60)
print("STARTING ALPHA SWEEP (linear VLM, no stall model)")
print("Range: -4 deg to 12 deg, step 2 deg")
print("=" * 60)

alphas_to_test = [float(a) for a in range(-4, 13, 2)]

results = {}
comp_geom_name = "VSPAEROComputeGeometry"

vsp.SetIntAnalysisInput(analysis_name, "StallModel", [0], 0)   # linear, no stall model

for alpha in alphas_to_test:
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [alpha], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [alpha], 0)

    vsp.SetAnalysisInputDefaults(comp_geom_name)
    vsp.ExecAnalysis(comp_geom_name)

    results_id = vsp.ExecAnalysis(analysis_name)
    if results_id:
        vsp.WriteResultsCSVFile(results_id, "Aircraft_Stability.csv")
        vsp.DeleteAllResults()
        cl, cd, cm = get_cl_cd_cm()
        results[alpha] = (cl, cd, cm)
        print(f"  alpha = {alpha:5.1f} deg  ->  CL = {cl:.4f}   CD = {cd:.5f}   Cm = {cm:.4f}")
    else:
        results[alpha] = (None, None, None)

# ======================== results table ========================
print("\n" + "=" * 64)
print(f"{'ALPHA (deg)':<12} | {'CL':<12} | {'CD':<12} | {'CM':<12} | {'L/D':<8}")
print("=" * 64)

for alpha in alphas_to_test:
    cl, cd, cm = results.get(alpha, (None, None, None))
    if cl is not None:
        ld = cl / cd if cd else float('nan')
        print(f"{alpha:<12.1f} | {cl:<12.4f} | {cd:<12.5f} | {cm:<12.4f} | {ld:<8.2f}")
    else:
        print(f"{alpha:<12.1f} | {'Error':<12} | {'Error':<12} | {'Error':<12} | {'Error':<8}")
print("=" * 64)
print("\nResults saved in 'log_vsp_trim.txt'")

# ======================== plots (one PNG each) ========================
import math
import matplotlib.pyplot as plt

a = alphas_to_test


def _v(x):
    return x if x is not None else math.nan


cl = [_v(results[x][0]) for x in a]
cd = [_v(results[x][1]) for x in a]
cm = [_v(results[x][2]) for x in a]
ld = [c / d if d else math.nan for c, d in zip(cl, cd)]

# ======================== curve fits -> derived aerodynamic parameters ========================
import numpy as np

aa = np.array(a, float)
CLa = np.array(cl, float)
CDa = np.array(cd, float)
CMa = np.array(cm, float)
deg2rad = np.pi / 180.0

# linear region (pre-stall, near the operating point) for the slopes
mask = aa <= 6.0

# lift curve: CL = CL_alpha*alpha + CL0
sCL, iCL = np.polyfit(aa[mask], CLa[mask], 1)
CL_alpha_deg = sCL
alpha_L0 = -iCL / sCL
CL0 = iCL

# pitching moment: Cm = Cm_alpha*alpha + Cm0
sCM, iCM = np.polyfit(aa[mask], CMa[mask], 1)
Cm_alpha_deg = sCM

# parabolic drag polar: CD = CD0 + k*CL^2
k_polar, CD0 = np.polyfit(CLa ** 2, CDa, 1)

# efficiency
LDmax = np.nanmax(ld)
a_LDmax = aa[int(np.nanargmax(ld))]

# static stability about the reference CG
SM = -Cm_alpha_deg / CL_alpha_deg   # static margin (fraction of MAC)
x_np = X_cg + SM * mean_chord        # neutral point

print("\n" + "=" * 56)
print("DERIVED AERODYNAMIC PARAMETERS (linear fit, alpha <= 6 deg)")
print("=" * 56)
print(f"  CL_alpha      = {CL_alpha_deg:.4f} /deg  ({CL_alpha_deg/deg2rad:.3f} /rad)")
print(f"  alpha_L0      = {alpha_L0:.2f} deg")
print(f"  CL0           = {CL0:.4f}")
print(f"  CD0           = {CD0:.5f}")
print(f"  k (CD=CD0+kCL^2)= {k_polar:.5f}")
print(f"  (L/D)max      = {LDmax:.2f}  at alpha = {a_LDmax:.1f} deg")
print(f"  Cm_alpha      = {Cm_alpha_deg:.4f} /deg  ({Cm_alpha_deg/deg2rad:.3f} /rad)")
print(f"  Static margin = {SM*100:.1f}% MAC   |   Neutral point = {x_np:.3f} m")
print("=" * 56)


def _save(fig, name):
    fig.savefig(name, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {name}")


# 1) Lift curve: CL vs alpha
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(a, cl, 'o-', color='navy')
ax.axhline(0.0, color='gray', lw=0.8)
ax.set_xlabel(r'$\alpha$ (deg)'); ax.set_ylabel(r'$C_L$')
ax.set_title('Lift curve'); ax.grid(True, ls='--', alpha=0.5)
_save(fig, 'lift_curve.png')

# 2) Drag polar: CL vs CD
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(cd, cl, 'o-', color='firebrick')
ax.axhline(0.0, color='gray', lw=0.8)
ax.set_xlabel(r'$C_D$'); ax.set_ylabel(r'$C_L$')
ax.set_title('Drag polar'); ax.grid(True, ls='--', alpha=0.5)
_save(fig, 'drag_polar.png')

# 3) Pitching moment: Cm vs alpha
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(a, cm, 'o-', color='teal')
ax.axhline(0.0, color='gray', lw=0.8)
ax.set_xlabel(r'$\alpha$ (deg)'); ax.set_ylabel(r'$C_m$')
ax.set_title('Pitching moment'); ax.grid(True, ls='--', alpha=0.5)
_save(fig, 'pitching_moment.png')

# 4) Aerodynamic efficiency: L/D vs alpha
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(a, ld, 'o-', color='darkorange')
ax.axhline(0.0, color='gray', lw=0.8)
ax.set_xlabel(r'$\alpha$ (deg)'); ax.set_ylabel(r'$L/D$')
ax.set_title('Aerodynamic efficiency'); ax.grid(True, ls='--', alpha=0.5)
_save(fig, 'lift_to_drag.png')

plt.show()
