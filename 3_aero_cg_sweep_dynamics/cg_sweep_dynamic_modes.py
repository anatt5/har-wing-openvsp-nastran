"""
Centre-of-gravity sweep: trim + dynamic modes (Etkin) vs x_cg.

Builds the aircraft geometry (main wing, horizontal tail, vertical tail,
control surfaces) at the fixed mesh density (d=50 panels/m) selected by the
separate VSPAERO mesh-convergence study, then for every x_cg in the sweep:
  1) rebuilds the point-mass model forcing the CG to x_cg (inertias
     recomputed), and saves the .vsp3 WITH the point masses
     (AC_cg_<xcg>.vsp3);
  2) trims the aircraft (alpha, elevator);
  3) runs the VSPAERO stability analysis at the trim point (AC.flt);
  4) assembles the Etkin longitudinal/lateral state-space matrices and
     extracts the dynamic modes (short period, phugoid, dutch roll, roll,
     spiral).
The short-period and phugoid modes are then tracked across the sweep using
the Modal Assurance Criterion (MAC) on their eigenvectors, since their
eigenvalues can cross in magnitude near the neutral point.

Outputs:
  AC.vsp3                    -- clean geometry (fixed d=50 mesh, no mass model)
  AC_cg_<xcg>.vsp3            -- geometry + point-mass model for each x_cg
  cg_sweep_results.csv        -- trim, static margin, inertias, modal results
  cg_sweep_log.txt            -- full console log
  cg_trim.png, cg_static_margin.png, cg_inertia.png,
  cg_longitudinal_modes.png, cg_stability_boundary.png,
  cg_lateral_modes.png, cg_root_locus_long.png, cg_root_locus_lat.png

Run:  python cg_sweep_dynamic_modes.py
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

# elevators: symmetric, deflected to a fixed trim angle
g_elevators_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Elevators_Group", g_elevators_index)
vsp.AddSelectedToCSGroup([5, 6], g_elevators_index)
set_surface_gain_by_index(2, 1.0, -1.0)

container_id = vsp.FindContainer("VSPAEROSettings", 0)
group_name = f"ControlSurfaceGroup_{g_elevators_index}"
parm_id = vsp.FindParm(container_id, "DeflectionAngle", group_name)
elevator_deflection = 2.7358
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
# Uniform mesh chosen as grid-independent: wing 25/81, HTP 29/20, VTP 37/12.
# Num W = Tess_W ("Shape"); Num U = SectTess_U (set on every xsec to be safe).
for _gid, _nw, _nu in ((wing_id, 25, 81), (htail_id, 29, 20), (vtp_id, 37, 12)):
    vsp.SetParmVal(_gid, "Tess_W", "Shape", _nw)
    _surf = vsp.GetXSecSurf(_gid, 0)
    for _i in range(vsp.GetNumXSec(_surf)):
        vsp.SetParmVal(vsp.GetXSecParm(vsp.GetXSec(_surf, _i), "SectTess_U"), float(_nu))
vsp.Update()

########################################################################
# 8) CENTRE-OF-GRAVITY SWEEP -> trim + dynamic modes for each x_cg
########################################################################
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ----- log everything to a file as well as the console -----
class _Logger(object):
    def __init__(self, filename="cg_sweep_log.txt"):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message); self.log.write(message); self.log.flush()

    def flush(self):
        self.terminal.flush(); self.log.flush()


sys.stdout = _Logger("cg_sweep_log.txt")

# ======================== reference constants ========================
m_total = 20.0    # total mass (kg)
g = 9.80665
rho = 1.225
S = 0.8           # Sref (m2)
bref = 3.2
mac = 0.25321     # cref / MAC (m)
v = 30.0          # Vinf (m/s) -- TU-Flex cruise
a_sound = 340.29
mu = 1.789e-5
Z_cg = 0.023      # CG z used in the VSPAERO runs

W = m_total * g
q_dyn = 0.5 * rho * v ** 2
CL_objective = W / (q_dyn * S)
mach = v / a_sound
re_cref = rho * v * mac / mu

g_elevators_index = 2   # 0 ailerons, 2 flaps, 3 elevators... (elevator group)
tolerance = 0.002
max_iter = 12

BASE_VSP = "AC.vsp3"    # clean geometry (no mass model) used as the base each iteration

# component masses (TU-Flex): fuselage 3.15, empennage 1.05, wings 3.98, rest = systems
M_SYS = m_total - (3.98 + 0.556 + 0.494 + 3.15)   # ~ 11.82 kg
# density reference (density 40 -> these surface masses, from AC_MassProps.txt)
DENS_REF = (("main wing", 0.64523, 3.98),
            ("HTP", 0.17188, 0.556),
            ("VTP", 0.15246, 0.494))

# save the clean base geometry (built above, with the d=50 mesh)
vsp.WriteVSPFile(BASE_VSP)
print(f"base geometry saved: {BASE_VSP}")


# ============================== helpers ==============================
def parse_massprops(path="AC_MassProps.txt"):
    d = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.endswith("Total Mass"):
                d["m"] = float(s.split()[0])
            elif s.endswith("Center of Gravity"):
                p = s.split(); d["cg"] = (float(p[0]), float(p[1]), float(p[2]))
            elif s.endswith("Ixx, Iyy, Izz"):
                p = s.split(); d["I1"] = (float(p[0]), float(p[1]), float(p[2]))
            elif s.endswith("Ixy, Ixz, Iyz"):
                p = s.split(); d["I2"] = (float(p[0]), float(p[1]), float(p[2]))
    return d


def build_mass_model(x_cg_ref):
    """Reload base geometry, build the point-mass model forcing CG = x_cg_ref,
    save the .vsp3 with point masses, return the MassProp dict."""
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(BASE_VSP)
    vsp.Update()

    # 1) lifting-surface masses via density
    for name, cur, tgt in DENS_REF:
        gid = vsp.FindGeomsWithName(name)[0]
        vsp.SetParmVal(gid, "Density", "Mass_Props", 40.0 * tgt / cur)
    vsp.Update()

    def add_pm(name, mass, x, y=0.0, z=0.0):
        bid = vsp.AddGeom("BLANK")
        vsp.SetGeomName(bid, name)
        vsp.SetParmVal(bid, "PointMass", "Mass_Props", mass)
        vsp.SetParmVal(bid, "X_Rel_Location", "XForm", x)
        vsp.SetParmVal(bid, "Y_Rel_Location", "XForm", y)
        vsp.SetParmVal(bid, "Z_Rel_Location", "XForm", z)
        vsp.Update()
        return bid

    # 2) fuselage (3 point masses) + 3) systems (forced to give CG = x_cg_ref)
    add_pm("fuselage_fwd", 1.05, 0.60)
    add_pm("fuselage_mid", 1.05, 1.50)
    add_pm("fuselage_aft", 1.05, 2.40)
    x0 = 1.05
    sid = add_pm("systems", M_SYS, x0)

    pm_set = set(vsp.FindGeoms())   # surfaces + point masses (mesh added by MassProp later)

    def _massprops():
        vsp.SetAnalysisInputDefaults("MassProp")
        try:
            vsp.SetIntAnalysisInput("MassProp", "NumMassSlices", [40])
        except Exception:
            pass
        vsp.ExecAnalysis("MassProp")
        return parse_massprops()

    mp = _massprops()
    dx = (x_cg_ref - mp["cg"][0]) * mp["m"] / M_SYS   # one-step CG correction
    vsp.SetParmVal(sid, "X_Rel_Location", "XForm", x0 + dx)
    vsp.Update()
    mp = _massprops()

    # remove the temporary MassProp mesh, keep the point masses, then save
    for gid in vsp.FindGeoms():
        if gid not in pm_set:
            try:
                vsp.DeleteGeom(gid)
            except Exception:
                pass
    vsp.Update()
    path = f"AC_cg_{x_cg_ref:.3f}.vsp3"
    vsp.WriteVSPFile(path)
    print(f"   mass model saved (with point masses): {path}")
    return mp


def setup_vspaero_inputs(x_cg):
    an = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(an)
    vsp.SetDoubleAnalysisInput(an, "Sref", [S], 0)
    vsp.SetDoubleAnalysisInput(an, "bref", [bref], 0)
    vsp.SetDoubleAnalysisInput(an, "cref", [mac], 0)
    vsp.SetDoubleAnalysisInput(an, "Xcg", [x_cg], 0)
    vsp.SetDoubleAnalysisInput(an, "Ycg", [0.0], 0)
    vsp.SetDoubleAnalysisInput(an, "Zcg", [Z_cg], 0)
    vsp.SetDoubleAnalysisInput(an, "Rho", [rho], 0)
    vsp.SetDoubleAnalysisInput(an, "Vinf", [v], 0)
    vsp.SetDoubleAnalysisInput(an, "MachStart", [mach], 0)
    vsp.SetDoubleAnalysisInput(an, "ReCref", [re_cref], 0)
    vsp.SetIntAnalysisInput(an, "Symmetry", [0], 0)
    vsp.SetIntAnalysisInput(an, "WakeNumIter", [5], 0)
    vsp.SetIntAnalysisInput(an, "NumWakeNodes", [64], 0)
    vsp.SetDoubleAnalysisInput(an, "WakeRelax", [0.5], 0)
    vsp.SetDoubleAnalysisInput(an, "CoreSizeFactor", [1.0], 0)
    vsp.SetIntAnalysisInput(an, "StallModel", [0], 0)
    vsp.SetDoubleAnalysisInput(an, "ThinGeomSet", [0], 0)
    vsp.SetDoubleAnalysisInput(an, "GeomSet", [0], 0)
    vsp.SetIntAnalysisInput(an, "AlphaNpts", [1], 0)
    vsp.SetIntAnalysisInput(an, "BetaNpts", [1], 0)
    vsp.SetDoubleAnalysisInput(an, "BetaStart", [0.0], 0)
    vsp.SetDoubleAnalysisInput(an, "BetaEnd", [0.0], 0)
    cid = vsp.FindContainer("VSPAEROSettings", 0)
    return vsp.FindParm(cid, "DeflectionAngle", f"ControlSurfaceGroup_{g_elevators_index}")


def _read_cl_cm(filename="Aircraft_Stability.csv"):
    cl = cm_ = None
    if not os.path.exists(filename):
        return None, None
    with open(filename) as f:
        for line in f:
            s = line.strip()
            if s.startswith("CLtot,"):
                p = s.split(",")
                if len(p) > 1:
                    cl = float(p[-1])
            elif s.startswith("CMytot,"):
                p = s.split(",")
                if len(p) > 1:
                    cm_ = float(p[-1])
    return cl, cm_


def execute_simulation(alpha, delta_e, parm_id):
    an = "VSPAEROSweep"
    vsp.SetDoubleAnalysisInput(an, "AlphaStart", [alpha], 0)
    vsp.SetDoubleAnalysisInput(an, "AlphaEnd", [alpha], 0)
    vsp.SetParmVal(parm_id, delta_e)
    vsp.Update()
    vsp.SetAnalysisInputDefaults("VSPAEROComputeGeometry")
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "ThinGeomSet", [0], 0)
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", [-1], 0)
    vsp.ExecAnalysis("VSPAEROComputeGeometry")
    rid = vsp.ExecAnalysis(an)
    if not rid:
        return None, None
    vsp.WriteResultsCSVFile(rid, "Aircraft_Stability.csv")
    vsp.DeleteAllResults()
    return _read_cl_cm()


def run_trim(parm_id):
    cl0, cm0 = execute_simulation(0.0, 0.0, parm_id)
    cl2, _ = execute_simulation(2.0, 0.0, parm_id)
    _, cmd5 = execute_simulation(0.0, 5.0, parm_id)
    CL_alpha = (cl2 - cl0) / 2.0
    Cm_delta = (cmd5 - cm0) / 5.0
    a, de, trimmed = 2.0, 0.0, False
    for _ in range(max_iter):
        CL, CM = execute_simulation(a, de, parm_id)
        if CL is None:
            break
        eCL = CL_objective - CL
        eCm = -CM
        if abs(eCL) < tolerance and abs(eCm) < tolerance:
            trimmed = True
            break
        a += eCL / CL_alpha
        de += eCm / Cm_delta
        a = max(-4.0, min(a, 10.0))
        de = max(-20.0, min(de, 20.0))
    return a, de, CL_alpha, Cm_delta, trimmed


def parse_flt_file(file_path, target_aoa):
    with open(file_path, "r") as f:
        content = f.read()
    for block in content.split("#" * 18):
        if not block.strip():
            continue
        data = {}
        for line in block.split("\n"):
            if ":" in line:
                k = line.split(":")[0].strip()
                try:
                    data[k] = float(line.split(":")[1].split()[0])
                except (ValueError, IndexError):
                    continue
        if "ALPHA_o" in data and np.isclose(data["ALPHA_o"], target_aoa, atol=0.1):
            return data
    raise ValueError(f"AoA = {target_aoa} deg not found in {file_path}")


def run_stability(a_trim, de_trim, parm_id):
    an = "VSPAEROSweep"
    vsp.SetParmVal(parm_id, de_trim)
    vsp.SetDoubleAnalysisInput(an, "AlphaStart", [a_trim], 0)
    vsp.SetDoubleAnalysisInput(an, "AlphaEnd", [a_trim], 0)
    vsp.SetIntAnalysisInput(an, "AlphaNpts", [1], 0)
    vsp.SetIntAnalysisInput(an, "UnsteadyType", [1], 0)   # stability -> writes AC.flt / AC.stab
    vsp.Update()
    vsp.SetAnalysisInputDefaults("VSPAEROComputeGeometry")
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "ThinGeomSet", [0], 0)
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", [-1], 0)
    vsp.ExecAnalysis("VSPAEROComputeGeometry")
    vsp.ExecAnalysis(an)
    return parse_flt_file("AC.flt", a_trim)


def build_state_matrices(flt, mp, a_trim):
    CL0 = flt.get("CLo", 0.0); CD0 = flt.get("CDo", 0.0)
    CLa = flt.get("CL_alpha", 0.0); CDa = flt.get("CD_alpha", 0.0); CLq = flt.get("CL_q", 0.0)
    Cma = flt.get("Cm_alpha", 0.0); Cmq = flt.get("Cm_q", 0.0)
    CLu = flt.get("CL_mach", 0.0); CDu = flt.get("CD_mach", 0.0); Cmu = flt.get("Cm_mach", 0.0)
    CYb = flt.get("CY_beta", 0.0); CYp = flt.get("CY_p", 0.0); CYr = flt.get("CY_r", 0.0)
    Clb = flt.get("Cl_beta", 0.0); Clp = flt.get("Cl_p", 0.0); Clr = flt.get("Cl_r", 0.0)
    Cnb = flt.get("Cn_beta", 0.0); Cnp = flt.get("Cn_p", 0.0); Cnr = flt.get("Cn_r", 0.0)

    m = mp["m"]; Ix = mp["I1"][0]; Iy = mp["I1"][1]; Iz = mp["I1"][2]; Ixz = mp["I2"][1]
    c = mac; b = bref; uo = v
    th = np.deg2rad(a_trim)

    Cxu = -(CDu + 2 * CD0); Cxa = CL0 - CDa
    Czu = -(CLu + 2 * CL0); Cza = -(CLa + CD0); Czq = -CLq

    Xu = 0.5 * rho * uo * S * Cxu; Xw = 0.5 * rho * uo * S * Cxa
    Zu = 0.5 * rho * uo * S * Czu; Zw = 0.5 * rho * uo * S * Cza; Zq = 0.25 * rho * uo * c * S * Czq
    Mu = 0.5 * rho * uo * c * S * Cmu
    Mw = 0.5 * rho * uo * c * S * Cma   # Etkin M_w (uo factor included)
    Mq = 0.25 * rho * uo * c * c * S * Cmq

    A_long = np.array([
        [Xu / m, Xw / m, 0.0, -g * np.cos(th)],
        [Zu / m, Zw / m, (Zq + m * uo) / m, -g * np.sin(th)],
        [Mu / Iy, Mw / Iy, Mq / Iy, 0.0],
        [0.0, 0.0, 1.0, 0.0]])

    Yv = 0.5 * rho * uo * S * CYb; Yp = 0.25 * rho * uo * S * b * CYp; Yr = 0.25 * rho * uo * S * b * CYr
    Lv = 0.5 * rho * uo * S * b * Clb; Lp = 0.25 * rho * uo * S * b * b * Clp; Lr = 0.25 * rho * uo * S * b * b * Clr
    Nv = 0.5 * rho * uo * S * b * Cnb; Np = 0.25 * rho * uo * S * b * b * Cnp; Nr = 0.25 * rho * uo * S * b * b * Cnr
    Ixp = (Ix * Iz - Ixz ** 2) / Iz; Izp = (Ix * Iz - Ixz ** 2) / Ix; Izxp = Ixz / (Ix * Iz - Ixz ** 2)

    A_lat = np.array([
        [Yv / m, Yp / m, (Yr / m) - uo, g * np.cos(th)],
        [Lv / Ixp + Izxp * Nv, Lp / Ixp + Izxp * Np, Lr / Ixp + Izxp * Nr, 0.0],
        [Izxp * Lv + Nv / Izp, Izxp * Lp + Np / Izp, Izxp * Lr + Nr / Izp, 0.0],
        [0.0, 1.0, np.tan(th), 0.0]])

    SM = -Cma / CLa if CLa != 0 else np.nan
    return A_long, A_lat, SM, CLa, Cma


def _osc(ev):
    """Return list of (wn, zeta, eig) for oscillatory eigenvalues (imag>0), sorted wn desc."""
    out = []
    for e in ev:
        if e.imag > 1e-6:
            wn = abs(e)
            out.append((wn, -e.real / wn, e))
    out.sort(key=lambda t: -t[0])
    return out


def extract_long(ev):
    osc = _osc(ev)
    if len(osc) >= 2:
        sp, ph = osc[0], osc[-1]
    elif len(osc) == 1:
        # Only one complex-conjugate pair left: the short period has already
        # split into two real roots (past the neutral point), so the
        # remaining oscillatory pair is the phugoid, not the short period.
        sp, ph = (np.nan, np.nan, None), osc[0]
    else:
        sp = ph = (np.nan, np.nan, None)
    return sp[0], sp[1], ph[0], ph[1], float(np.max(ev.real))


def extract_lat(ev):
    osc = _osc(ev)
    dr = osc[0] if osc else (np.nan, np.nan, None)
    reals = sorted([e.real for e in ev if abs(e.imag) <= 1e-6])
    roll = reals[0] if len(reals) >= 1 else np.nan
    spiral = reals[-1] if len(reals) >= 2 else np.nan
    return dr[0], dr[1], roll, spiral


def _mac(a, b):
    """Modal Assurance Criterion between two eigenvectors: 1 = same mode
    shape, 0 = unrelated shapes. Scale/phase-independent."""
    num = abs(np.vdot(a, b)) ** 2
    den = np.vdot(a, a).real * np.vdot(b, b).real
    return num / den


def _best_pair_score(vecs_a, vecs_b):
    """Best total MAC when matching two candidate eigenvectors against two
    reference eigenvectors, trying both possible pairings."""
    m11 = _mac(vecs_a[0], vecs_b[0]); m22 = _mac(vecs_a[1], vecs_b[1])
    m12 = _mac(vecs_a[0], vecs_b[1]); m21 = _mac(vecs_a[1], vecs_b[0])
    return max(m11 + m22, m12 + m21)


def _group_blocks(ev):
    """Split eigenvalues into indivisible blocks: a complex-conjugate pair
    is one block of 2 (its two eigenvalues are one physical mode and must
    never end up in different families); a real eigenvalue is its own
    block of 1."""
    ev = list(ev)
    used = [False] * len(ev)
    blocks = []
    for i, e in enumerate(ev):
        if used[i]:
            continue
        if abs(e.imag) <= 1e-6:
            blocks.append([i]); used[i] = True
        else:
            for j in range(len(ev)):
                if not used[j] and j != i and abs(ev[j] - np.conj(e)) < 1e-6:
                    blocks.append([i, j]); used[i] = used[j] = True
                    break
            else:
                blocks.append([i]); used[i] = True
    return blocks


def track_longitudinal_modes(results):
    """Track the short-period and phugoid modes across the CG sweep using
    the Modal Assurance Criterion (MAC) on the eigenvectors, not just the
    eigenvalues' position in the complex plane.

    Why: near the neutral point the two modes' eigenvalues can become close
    in magnitude, so distance-based tracking becomes ambiguous (and can even
    split a genuine complex-conjugate pair between the two families). But
    the short period (mostly w/q/theta motion, u roughly constant) and the
    phugoid (mostly u/theta motion, alpha roughly constant) have very
    different *mode shapes* -- a difference that survives the eigenvalue
    crossover. A complex-conjugate pair is always kept together as a single
    block, since it represents one physical mode.

    Returns sp_wn, sp_z, ph_wn, ph_z (NaN once that mode's tracked
    eigenvalue is real) and, for each x_cg, the raw eigenvalues tracked as
    short period / phugoid (for the root-locus plot).
    """
    n = len(results)
    sp_wn = np.full(n, np.nan); sp_z = np.full(n, np.nan)
    ph_wn = np.full(n, np.nan); ph_z = np.full(n, np.nan)
    sp_ev_list = [None] * n; ph_ev_list = [None] * n

    # Seed at the first x_cg, where SP (high frequency) and phugoid (near
    # origin) are unambiguous by magnitude alone.
    ev0, evec0 = results[0]["ev_long"], results[0]["evec_long"]
    order0 = list(np.argsort(-np.abs(ev0)))
    prev_sp_idx, prev_ph_idx = order0[:2], order0[2:]
    prev_sp_vecs = [evec0[:, k] for k in prev_sp_idx]
    prev_ph_vecs = [evec0[:, k] for k in prev_ph_idx]

    for i, r in enumerate(results):
        ev, evec = r["ev_long"], r["evec_long"]
        blocks = _group_blocks(ev)

        # Try every way to split the blocks into two groups of exactly 2
        # eigenvalues each (never breaking a block apart), and keep the
        # split+assignment that maximises total MAC similarity to the
        # previous step's SP/PH mode shapes.
        best = None
        for r_bits in range(1 << len(blocks)):
            group_a = [blocks[b] for b in range(len(blocks)) if (r_bits >> b) & 1]
            group_b = [blocks[b] for b in range(len(blocks)) if not (r_bits >> b) & 1]
            flat_a = [j for g in group_a for j in g]
            flat_b = [j for g in group_b for j in g]
            if len(flat_a) != 2 or len(flat_b) != 2:
                continue
            vecs_a = [evec[:, j] for j in flat_a]
            vecs_b = [evec[:, j] for j in flat_b]
            for sp_idx, sp_vecs, ph_idx, ph_vecs in (
                (flat_a, vecs_a, flat_b, vecs_b),
                (flat_b, vecs_b, flat_a, vecs_a),
            ):
                score = (_best_pair_score(sp_vecs, prev_sp_vecs)
                         + _best_pair_score(ph_vecs, prev_ph_vecs))
                if best is None or score > best[0]:
                    best = (score, sp_idx, ph_idx)

        _, sp_idx, ph_idx = best
        sp_e, ph_e = ev[sp_idx], ev[ph_idx]
        sp_ev_list[i] = sp_e
        ph_ev_list[i] = ph_e
        for e in sp_e:
            if e.imag > 1e-6:
                wn = abs(e); sp_wn[i] = wn; sp_z[i] = -e.real / wn
        for e in ph_e:
            if e.imag > 1e-6:
                wn = abs(e); ph_wn[i] = wn; ph_z[i] = -e.real / wn

        prev_sp_vecs = [evec[:, j] for j in sp_idx]
        prev_ph_vecs = [evec[:, j] for j in ph_idx]

    return sp_wn, sp_z, ph_wn, ph_z, sp_ev_list, ph_ev_list


# ============================== CG SWEEP ==============================
xcgs = np.round(np.arange(1.20, 1.40 + 1e-9, 0.025), 3)
results = []

for xcg in xcgs:
    print("\n" + "#" * 70)
    print(f"#  x_cg = {xcg:.3f} m")
    print("#" * 70)

    mp = build_mass_model(xcg)
    print(f"   m={mp['m']:.3f} kg  CG_x={mp['cg'][0]:.4f}  "
          f"Ixx={mp['I1'][0]:.3f} Iyy={mp['I1'][1]:.3f} Izz={mp['I1'][2]:.3f} Ixz={mp['I2'][1]:.3f}")

    # reload pristine geometry for the aerodynamic / stability runs
    vsp.ClearVSPModel(); vsp.ReadVSPFile(BASE_VSP); vsp.Update()
    parm_id = setup_vspaero_inputs(xcg)

    a_trim, de_trim, CLa_fit, Cmd_fit, trimmed = run_trim(parm_id)
    if not trimmed:
        print(f"   [!] trim not converged at x_cg={xcg:.3f} -> skipped")
        continue
    print(f"   TRIM: alpha={a_trim:.4f} deg | elevator={de_trim:.4f} deg")

    flt = run_stability(a_trim, de_trim, parm_id)
    A_long, A_lat, SM, CLa, Cma = build_state_matrices(flt, mp, a_trim)
    ev_long, evec_long = np.linalg.eig(A_long)
    ev_lat = np.linalg.eigvals(A_lat)

    sp_wn, sp_z, ph_wn, ph_z, long_maxreal = extract_long(ev_long)
    dr_wn, dr_z, roll, spiral = extract_lat(ev_lat)

    print(f"   SM={SM*100:6.2f}%MAC | SP: {sp_wn:.3f} rad/s z={sp_z:.3f} "
          f"({sp_wn/(2*np.pi):.3f} Hz) | Phug: {ph_wn:.3f} rad/s z={ph_z:.3f}")

    results.append(dict(xcg=xcg, alpha=a_trim, dele=de_trim, SM=SM,
                         Ixx=mp["I1"][0], Iyy=mp["I1"][1], Izz=mp["I1"][2], Ixz=mp["I2"][1],
                         ev_long=ev_long, evec_long=evec_long, ev_lat=ev_lat,
                         sp_wn=sp_wn, sp_z=sp_z, ph_wn=ph_wn, ph_z=ph_z, long_maxreal=long_maxreal,
                         dr_wn=dr_wn, dr_z=dr_z, roll=roll, spiral=spiral))

# Re-classify short period vs phugoid by continuity across the sweep (see
# track_longitudinal_modes docstring) instead of the per-point heuristic
# used inline in the loop above, which is unreliable near the neutral point.
sp_wn_trk, sp_z_trk, ph_wn_trk, ph_z_trk, sp_ev_list, ph_ev_list = track_longitudinal_modes(results)

# Once the classical two-mode picture collapses completely (no
# complex-conjugate pairs left at all -- both SP and phugoid momentarily
# real, as happens right at the neutral point), any oscillatory pair that
# reappears further aft comes from real roots re-merging (e.g. the
# divergent short-period root interacting with a phugoid root) and cannot
# be cleanly attributed to either historical mode. Freeze both curves from
# that point on instead of resurrecting a label.
collapsed = False
for i, r in enumerate(results):
    if len(_osc(r["ev_long"])) == 0:
        collapsed = True
    r["collapsed"] = collapsed
    if collapsed:
        sp_wn_trk[i] = sp_z_trk[i] = np.nan
        ph_wn_trk[i] = ph_z_trk[i] = np.nan

for i, r in enumerate(results):
    r["sp_wn"], r["sp_z"] = sp_wn_trk[i], sp_z_trk[i]
    r["ph_wn"], r["ph_z"] = ph_wn_trk[i], ph_z_trk[i]
    r["sp_ev"], r["ph_ev"] = sp_ev_list[i], ph_ev_list[i]

# ============================== SAVE TABLE ==============================
arr = {k: np.array([r[k] for r in results]) for k in
       ("xcg", "alpha", "dele", "SM", "Iyy", "sp_wn", "sp_z", "ph_wn", "ph_z",
        "long_maxreal", "dr_wn", "dr_z", "roll", "spiral")}

with open("cg_sweep_results.csv", "w") as f:
    f.write("xcg,alpha_trim,elevator_trim,SM_%MAC,Iyy,SP_freq_Hz,SP_zeta,"
             "Phug_freq_Hz,Phug_zeta,long_max_real,DR_freq_Hz,DR_zeta,roll_eig,spiral_eig\n")
    for r in results:
        f.write(f"{r['xcg']:.3f},{r['alpha']:.4f},{r['dele']:.4f},{r['SM']*100:.3f},{r['Iyy']:.4f},"
                f"{r['sp_wn']/(2*np.pi):.4f},{r['sp_z']:.4f},"
                f"{r['ph_wn']/(2*np.pi):.4f},{r['ph_z']:.4f},{r['long_maxreal']:.4f},"
                f"{r['dr_wn']/(2*np.pi):.4f},{r['dr_z']:.4f},{r['roll']:.4f},{r['spiral']:.4f}\n")
print("\nresults table saved: cg_sweep_results.csv")

# ============================== PLOTS ==============================
NP_xcg = None   # where SM crosses zero (neutral point), for reference lines
if np.any(arr["SM"] > 0) and np.any(arr["SM"] < 0):
    NP_xcg = np.interp(0.0, arr["SM"][::-1] * 100, arr["xcg"][::-1])


def _np_line(ax):
    if NP_xcg is not None:
        ax.axvline(NP_xcg, color="r", ls=":", lw=1.2, label=f"NP $\\approx$ {NP_xcg:.3f} m")


# 1) Trim alpha and elevator vs CG
fig, ax1 = plt.subplots(figsize=(7, 4.5))
ax1.plot(arr["xcg"], arr["alpha"], "o-", color="tab:blue", label=r"$\alpha_{trim}$")
ax1.set_xlabel(r"$x_{cg}$ [m]"); ax1.set_ylabel(r"$\alpha_{trim}$ [deg]", color="tab:blue")
ax1.tick_params(axis="y", labelcolor="tab:blue")
ax2 = ax1.twinx()
ax2.plot(arr["xcg"], arr["dele"], "s--", color="tab:red", label=r"$\delta_{e,trim}$")
ax2.set_ylabel(r"$\delta_{e,trim}$ [deg]", color="tab:red")
ax2.tick_params(axis="y", labelcolor="tab:red")
_np_line(ax1)
ax1.set_title("Trim condition vs centre of gravity"); ax1.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("cg_trim.png", dpi=150); plt.close(fig)

# 2) Static margin vs CG
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(arr["xcg"], arr["SM"] * 100, "o-", color="tab:green")
ax.axhline(0, color="k", lw=0.8, ls="--")
_np_line(ax)
ax.set_xlabel(r"$x_{cg}$ [m]"); ax.set_ylabel("Static margin [\\% MAC]")
ax.set_title("Static margin vs centre of gravity"); ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("cg_static_margin.png", dpi=150); plt.close(fig)

# 3) Moments of inertia vs CG
fig, ax = plt.subplots(figsize=(7, 4.5))
for k, lab in (("Ixx", r"$I_{xx}$"), ("Iyy", r"$I_{yy}$"), ("Izz", r"$I_{zz}$")):
    ax.plot([r["xcg"] for r in results], [r[k] for r in results], "o-", label=lab)
ax.set_xlabel(r"$x_{cg}$ [m]"); ax.set_ylabel(r"Moment of inertia [kg$\cdot$m$^2$]")
ax.set_title("Moments of inertia vs centre of gravity"); ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("cg_inertia.png", dpi=150); plt.close(fig)

# 4) Longitudinal modes (short period + phugoid): frequency and damping vs CG
fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
axs[0, 0].plot(arr["xcg"], arr["sp_wn"] / (2 * np.pi), "o-", color="tab:blue")
axs[0, 0].set_ylabel("freq [Hz]"); axs[0, 0].set_title("Short period")
axs[1, 0].plot(arr["xcg"], arr["sp_z"], "o-", color="tab:blue")
axs[1, 0].set_ylabel(r"damping $\zeta$"); axs[1, 0].set_xlabel(r"$x_{cg}$ [m]")
axs[0, 1].plot(arr["xcg"], arr["ph_wn"] / (2 * np.pi), "s-", color="tab:orange")
axs[0, 1].set_ylabel("freq [Hz]"); axs[0, 1].set_title("Phugoid")
axs[1, 1].plot(arr["xcg"], arr["ph_z"], "s-", color="tab:orange")
axs[1, 1].set_ylabel(r"damping $\zeta$"); axs[1, 1].set_xlabel(r"$x_{cg}$ [m]")
for a_ in axs.flat:
    a_.grid(alpha=0.3); _np_line(a_)
fig.suptitle("Longitudinal modes vs centre of gravity")
fig.tight_layout(); fig.savefig("cg_longitudinal_modes.png", dpi=150); plt.close(fig)

# 5) Stability boundary: largest real part of the longitudinal eigenvalues vs CG
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(arr["xcg"], arr["long_maxreal"], "o-", color="tab:purple")
ax.axhline(0, color="k", lw=0.8, ls="--")
_np_line(ax)
ax.set_xlabel(r"$x_{cg}$ [m]"); ax.set_ylabel(r"max Re($\lambda$) longitudinal [1/s]")
ax.set_title("Longitudinal stability boundary vs centre of gravity")
ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("cg_stability_boundary.png", dpi=150); plt.close(fig)

# 6) Lateral-directional modes vs CG
fig, axs = plt.subplots(1, 3, figsize=(13, 4))
axs[0].plot(arr["xcg"], arr["dr_wn"] / (2 * np.pi), "o-", color="tab:blue", label="freq [Hz]")
axs[0].plot(arr["xcg"], arr["dr_z"], "s--", color="tab:red", label=r"$\zeta$")
axs[0].set_title("Dutch roll"); axs[0].set_xlabel(r"$x_{cg}$ [m]"); axs[0].legend()
axs[1].plot(arr["xcg"], arr["roll"], "o-", color="tab:green")
axs[1].set_title("Roll subsidence (eigenvalue)"); axs[1].set_xlabel(r"$x_{cg}$ [m]")
axs[1].set_ylabel("Re($\\lambda$) [1/s]")
axs[2].plot(arr["xcg"], arr["spiral"], "o-", color="tab:purple")
axs[2].axhline(0, color="k", lw=0.8, ls="--")
axs[2].set_title("Spiral mode (eigenvalue)"); axs[2].set_xlabel(r"$x_{cg}$ [m]")
axs[2].set_ylabel("Re($\\lambda$) [1/s]")
for a_ in axs:
    a_.grid(alpha=0.3)
fig.suptitle("Lateral-directional modes vs centre of gravity")
fig.tight_layout(); fig.savefig("cg_lateral_modes.png", dpi=150); plt.close(fig)

# 7) Longitudinal root locus (two colour bands: short period vs phugoid)
fig, ax = plt.subplots(figsize=(7.5, 6))
norm = plt.Normalize(arr["xcg"].min(), arr["xcg"].max())
cmap_sp, cmap_ph = cm.Blues, cm.Oranges
for r in results:
    # Use the same continuity-tracked short-period / phugoid labels as the
    # modes-vs-cg figure (see track_longitudinal_modes), so the two plots
    # agree and mode identity survives the crossover near the neutral point.
    shade = 0.35 + 0.55 * norm(r["xcg"])
    ax.scatter(r["sp_ev"].real, r["sp_ev"].imag, color=cmap_sp(shade), s=45, edgecolors="k", linewidths=0.4)
    ax.scatter(r["ph_ev"].real, r["ph_ev"].imag, color=cmap_ph(shade), s=45, edgecolors="k", linewidths=0.4)
ax.axvline(0, color="k", lw=0.9, ls="--")
ax.set_xlabel(r"Re($\lambda$) [1/s]"); ax.set_ylabel(r"Im($\lambda$) [1/s]")
ax.set_title("Longitudinal root locus vs $x_{cg}$"); ax.grid(alpha=0.3)
sm_sp = cm.ScalarMappable(norm=norm, cmap=cmap_sp); sm_sp.set_array([])
sm_ph = cm.ScalarMappable(norm=norm, cmap=cmap_ph); sm_ph.set_array([])
fig.colorbar(sm_sp, ax=ax, label=r"$x_{cg}$ [m] (short period)")
fig.colorbar(sm_ph, ax=ax, label=r"$x_{cg}$ [m] (phugoid)")
fig.tight_layout(); fig.savefig("cg_root_locus_long.png", dpi=150); plt.close(fig)

# 8) Lateral-directional root locus
fig, ax = plt.subplots(figsize=(7.5, 6))
for r in results:
    col = cm.viridis(norm(r["xcg"]))
    ax.scatter(r["ev_lat"].real, r["ev_lat"].imag, color=col, s=45, edgecolors="k", linewidths=0.4)
ax.axvline(0, color="k", lw=0.9, ls="--")
ax.set_xlabel(r"Re($\lambda$) [1/s]"); ax.set_ylabel(r"Im($\lambda$) [1/s]")
ax.set_title("Lateral-directional root locus vs $x_{cg}$"); ax.grid(alpha=0.3)
sm = cm.ScalarMappable(norm=norm, cmap="viridis"); sm.set_array([])
fig.colorbar(sm, ax=ax, label=r"$x_{cg}$ [m]")
fig.tight_layout(); fig.savefig("cg_root_locus_lat.png", dpi=150); plt.close(fig)

print("\n" + "=" * 60)
print("CG SWEEP COMPLETE")
print(f"  cases run: {len(results)}")
if NP_xcg is not None:
    print(f"  neutral point (SM=0) at x_cg = {NP_xcg:.3f} m")
print("  plots: cg_trim.png, cg_static_margin.png, cg_inertia.png,")
print("         cg_longitudinal_modes.png, cg_stability_boundary.png,")
print("         cg_lateral_modes.png, cg_root_locus_long.png, cg_root_locus_lat.png")
print("  data : cg_sweep_results.csv  |  models: AC_cg_<xcg>.vsp3")
print("=" * 60)
