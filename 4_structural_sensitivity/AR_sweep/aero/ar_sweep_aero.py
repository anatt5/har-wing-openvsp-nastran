"""
ASPECT RATIO SWEEP -- aerodynamic / flight-dynamics side.

Adapted from ../../../3_aero_cg_sweep_dynamics/cg_sweep_dynamic_modes.py (the CG sweep): all the
trim / VSPAERO-stability / Etkin state-matrix / mode-tracking machinery is
reused unchanged. What's different here:
  - The wing is built with the AREA + AR + TAPER driver group (instead of
    SPAN + AR + TAPER), with Area held fixed at the baseline wing area
    (0.8 m^2 total -> 0.4 m^2 per panel). This means increasing AR grows
    the span and shrinks the chord, holding wing area constant -- the
    classical "higher aspect ratio wing" mechanism that the TU-Flex papers
    cite for why HAR wings show lower aeroelastic-mode frequencies (a
    longer moment arm, not a thinner chord).
  - x_cg is held FIXED at 1.325 m (10.0% MAC static margin at the baseline
    AR, from cg_sweep_results.csv) for every AR point, to isolate the AR
    effect from the CG effect.
  - Sref/bref/cref are re-measured from the actual wing geometry at each
    AR (they are NOT constant here, unlike the CG sweep, since span/chord
    change with AR even though area doesn't).
  - The whole geometry build is wrapped in build_geometry(AR_mw) and
    called fresh for every AR value, instead of being built once before
    the loop.

Run:  python ar_sweep_aero.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import openvsp as vsp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
AIRFOIL_DIR = os.path.join(REPO_ROOT, "airfoils")

sys.path.insert(0, REPO_ROOT)
from utils.airfoil_utils import assign_airfoil_to_component

# ----- log everything to a file as well as the console -----
class _Logger(object):
    def __init__(self, filename="ar_sweep_log.txt"):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")
    def write(self, message):
        self.terminal.write(message); self.log.write(message); self.log.flush()
    def flush(self):
        self.terminal.flush(); self.log.flush()
sys.stdout = _Logger(os.path.join(HERE, "ar_sweep_log.txt"))

# ============================== reference constants ==============================
m_total = 20.0           # total mass (kg)
g       = 9.80665
rho     = 1.225
v       = 30.0            # Vinf (m/s) - TU-Flex cruise
a_sound = 340.29
mu      = 1.789e-5
Z_cg    = 0.023
X_CG_FIXED = 1.325         # 10.0% MAC static margin at baseline AR (cg_sweep_results.csv) -- held fixed across the AR sweep
S_FIXED = 0.8              # wing area held fixed (m^2) -- only span/chord change with AR
TAPER_MW = 0.66667
SWEEP_ANGLE = 28.0

W       = m_total * g
mach    = v / a_sound

g_elevators_index = 2      # 0 ailerons, 2 flaps, 3 elevators... (elevator group)
tolerance = 0.002
max_iter  = 12   # safety cap -- with warm-starting (below) it should converge in far fewer steps

BASE_VSP = os.path.join(HERE, "AC.vsp3")

# component masses (TU-Flex): fuselage 3.15, empennage 1.05, wings 3.98, rest = systems
M_SYS = m_total - (3.98 + 0.556 + 0.494 + 3.15)   # ~ 11.82 kg
# density reference (density 40 -> these surface masses, from AC_MassProps.txt)
# valid across the AR sweep since wing AREA is held fixed (only span/chord redistribute)
DENS_REF = (("main wing", 0.64523, 3.98),
            ("HTP",       0.17188, 0.556),
            ("VTP",       0.15246, 0.494))


# ============================== full geometry build ==============================
def build_geometry(AR_mw):
    """Build main wing (at the given panel AR, area held fixed) + HTP + VTP +
    control surfaces, exactly as in vspaerotrim_bucle.py, and save it as
    BASE_VSP. Returns (wing_id, htail_id, vtp_id).

    No fuselage surface: matches the original dynamic-modes methodology,
    which used only wing + HTP + VTP aerodynamically and represented the
    fuselage purely through the point-mass system (fuselage_fwd/mid/aft +
    systems, added in build_mass_model)."""
    vsp.VSPRenew()

    # ---- main wing: AREA + AR + TAPER driver, area fixed, AR swept ----
    wing_id = vsp.AddGeom("WING", "")
    vsp.SetGeomName(wing_id, "main wing")
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", 25)

    vsp.SetDriverGroup(wing_id, 1,
                    vsp.AREA_WSECT_DRIVER,
                    vsp.AR_WSECT_DRIVER,
                    vsp.TAPER_WSECT_DRIVER)
    vsp.Update()

    vsp.SetParmVal(vsp.GetParm(wing_id, "Area", "XSec_1"), S_FIXED / 2.0)  # per-panel area
    vsp.SetParmVal(vsp.GetParm(wing_id, "Aspect", "XSec_1"), AR_mw)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Taper", "XSec_1"), TAPER_MW)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Sweep", "XSec_1"), SWEEP_ANGLE)
    vsp.Update()

    vsp.SetParmVal(vsp.GetParm(wing_id, "X_Rel_Location", "XForm"), 0.7)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Location", "XForm"), 0)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Z_Rel_Location", "XForm"), 0)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Rotation", "XForm"), 0)
    vsp.Update()

    assign_airfoil_to_component(wing_id, "clark Y.dat", AIRFOIL_DIR)
    xsec_surf_id_mwing = vsp.GetXSecSurf(wing_id, 0)
    for i in (0, 1):
        xsec_id = vsp.GetXSec(xsec_surf_id_mwing, i)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "SectTess_U"), 81.0)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "Dihedral"), 6.0)
    vsp.SetParmVal(wing_id, "Density", "Mass_Props", 40)
    vsp.Update()

    # ---- HTP ----
    htail_id = vsp.AddGeom("WING")
    vsp.SetGeomName(htail_id, "HTP")
    vsp.SetParmVal(htail_id, "Tess_W", "Shape", 29)

    vsp.SetDriverGroup(htail_id, 1,
                    vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
    vsp.Update()
    vsp.SetParmVal(vsp.GetParm(htail_id, "Span",  "XSec_1"), 0.38)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Taper", "XSec_1"), 0.57143)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Root_Chord", "XSec_1"), 0.35)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Sweep", "XSec_1"), SWEEP_ANGLE)

    vsp.SetParmVal(vsp.GetParm(htail_id, "X_Rel_Location", "XForm"), 2.55)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Z_Rel_Location", "XForm"), 0.615)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Y_Rel_Location", "XForm"), 0.0)
    vsp.SetParmVal(vsp.GetParm(htail_id, "Y_Rel_Rotation", "XForm"), -3)
    vsp.Update()

    assign_airfoil_to_component(htail_id, "naca 0012.dat", AIRFOIL_DIR)
    xsec_surf_id_htp = vsp.GetXSecSurf(htail_id, 0)
    for i in (0, 1):
        xsec_id = vsp.GetXSec(xsec_surf_id_htp, i)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "SectTess_U"), 24.0)
    vsp.SetParmVal(htail_id, "Density", "Mass_Props", 40)
    xsec_1 = vsp.GetXSec(xsec_surf_id_htp, 1)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_1, "ThickChord"), 0.08)

    # ---- VTP (two sections) ----
    vtp_id = vsp.AddGeom("WING")
    vsp.SetGeomName(vtp_id, "VTP")
    vsp.SetParmVal(vtp_id, "Tess_W", "Shape", 37)
    vsp.InsertXSec(vtp_id, 1, vsp.XS_FOUR_SERIES)
    vsp.Update()

    vsp.SetDriverGroup(vtp_id, 1, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
    vsp.Update()
    vsp.SetParmVal(vtp_id, "Span", "XSec_1", 0.2)
    vsp.SetParmVal(vtp_id, "Taper", "XSec_1", 0.68966)
    vsp.SetParmVal(vtp_id, "Root_Chord", "XSec_1", 0.5)
    vsp.SetParmVal(vtp_id, "Sweep", "XSec_1", 50.0)

    vsp.SetDriverGroup(vtp_id, 2, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
    vsp.Update()
    vsp.SetParmVal(vtp_id, "Span", "XSec_2", 0.23)
    vsp.SetParmVal(vtp_id, "Taper", "XSec_2", 0.81034)
    vsp.SetParmVal(vtp_id, "Root_Chord", "XSec_2", 0.34483)
    vsp.SetParmVal(vtp_id, "Sweep", "XSec_2", 35.0)

    vsp.SetParmVal(vtp_id, "X_Rel_Rotation", "XForm", 90)
    vsp.SetParmVal(vtp_id, "Sym_Planar_Flag", "Sym", 0)
    vsp.SetParmVal(vtp_id, "X_Rel_Location", "XForm", 2.15)
    vsp.SetParmVal(vtp_id, "Z_Rel_Location", "XForm", 0.17)
    vsp.Update()

    assign_airfoil_to_component(vtp_id, "naca 0008.dat", AIRFOIL_DIR)
    xsec_surf_id = vsp.GetXSecSurf(vtp_id, 0)
    xsec_0 = vsp.GetXSec(xsec_surf_id, 0)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_0, "ThickChord"), 0.08)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_0, "SectTess_U"), 26.0)
    xsec_1v = vsp.GetXSec(xsec_surf_id, 1)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_1v, "ThickChord"), 0.12)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_1v, "SectTess_U"), 26.0)
    xsec_2v = vsp.GetXSec(xsec_surf_id, 2)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_2v, "ThickChord"), 0.08)
    vsp.SetParmVal(vsp.GetXSecParm(xsec_2v, "SectTess_U"), 26.0)
    vsp.Update()
    vsp.SetParmVal(vtp_id, "Density", "Mass_Props", 40)

    # ---- close geometry ----
    for gid in (wing_id, htail_id, vtp_id):
        vsp.SetParmVal(gid, "CapUMinOption", "EndCap", 1.0)
        vsp.SetParmVal(gid, "CapUMaxOption", "EndCap", 1.0)

    # ---- control surfaces ----
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
        all_settings = vsp.FindContainer("VSPAEROSettings", 0)
        all_parms = vsp.FindContainerParmIDs(all_settings)
        gain_pids = [pid for pid in all_parms if "Gain" in vsp.GetParmName(pid)]
        idx0, idx1 = group_index * 2, group_index * 2 + 1
        vsp.SetParmVal(gain_pids[idx0], target_val_0)
        vsp.SetParmVal(gain_pids[idx1], target_val_1)
        vsp.Update()

    g_ailerons_index = vsp.CreateVSPAEROControlSurfaceGroup()
    vsp.SetVSPAEROControlGroupName("Ailerons_Group", g_ailerons_index)
    vsp.AddSelectedToCSGroup([1, 2], g_ailerons_index)
    set_surface_gain_by_index(0, 1, 1)

    g_flaps_index = vsp.CreateVSPAEROControlSurfaceGroup()
    vsp.SetVSPAEROControlGroupName("Flaps_Group", g_flaps_index)
    vsp.AddSelectedToCSGroup([3, 4], g_flaps_index)
    set_surface_gain_by_index(1, 1.0, -1.0)

    g_elev_idx = vsp.CreateVSPAEROControlSurfaceGroup()
    vsp.SetVSPAEROControlGroupName("Elevators_Group", g_elev_idx)
    vsp.AddSelectedToCSGroup([5, 6], g_elev_idx)
    set_surface_gain_by_index(2, 1.0, -1.0)

    g_rudder_index = vsp.CreateVSPAEROControlSurfaceGroup()
    vsp.SetVSPAEROControlGroupName("Rudder_Group", g_rudder_index)
    vsp.AddSelectedToCSGroup([7], g_rudder_index)
    vsp.Update()

    # ---- mesh density (from the mesh convergence study) ----
    for _gid, _nw, _nu in ((wing_id, 25, 81), (htail_id, 29, 20), (vtp_id, 37, 12)):
        vsp.SetParmVal(_gid, "Tess_W", "Shape", _nw)
        _surf = vsp.GetXSecSurf(_gid, 0)
        for _i in range(vsp.GetNumXSec(_surf)):
            vsp.SetParmVal(vsp.GetXSecParm(vsp.GetXSec(_surf, _i), "SectTess_U"), float(_nu))
    vsp.Update()

    vsp.WriteVSPFile(BASE_VSP)
    return wing_id, htail_id, vtp_id


# ============================== mass / trim / stability (unchanged from CG sweep) ==============================
def parse_massprops(path):
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


def build_mass_model(x_cg_ref, AR_std):
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(BASE_VSP)
    vsp.Update()

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

    add_pm("fuselage_fwd", 1.05, 0.60)
    add_pm("fuselage_mid", 1.05, 1.50)
    add_pm("fuselage_aft", 1.05, 2.40)
    x0 = 1.05
    sid = add_pm("systems", M_SYS, x0)

    pm_set = set(vsp.FindGeoms())

    def _massprops():
        vsp.SetAnalysisInputDefaults("MassProp")
        try:
            vsp.SetIntAnalysisInput("MassProp", "NumMassSlices", [40])
        except Exception:
            pass
        vsp.ExecAnalysis("MassProp")
        return parse_massprops(os.path.join(HERE, "AC_MassProps.txt")
                                if os.path.exists(os.path.join(HERE, "AC_MassProps.txt"))
                                else "AC_MassProps.txt")

    mp = _massprops()
    dx = (x_cg_ref - mp["cg"][0]) * mp["m"] / M_SYS
    vsp.SetParmVal(sid, "X_Rel_Location", "XForm", x0 + dx)
    vsp.Update()
    mp = _massprops()

    for gid in vsp.FindGeoms():
        if gid not in pm_set:
            try:
                vsp.DeleteGeom(gid)
            except Exception:
                pass
    vsp.Update()
    path = os.path.join(HERE, f"AC_AR_{AR_std:.2f}.vsp3")
    vsp.WriteVSPFile(path)
    print(f"   mass model saved (with point masses): {path}")
    return mp


def setup_vspaero_inputs(x_cg, S, bref, mac, mach, re_cref):
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


def _read_cl_cm(filename):
    cl = cm_ = None
    if not os.path.exists(filename):
        return None, None
    with open(filename) as f:
        for line in f:
            s = line.strip()
            if s.startswith("CLtot,"):
                p = s.split(",")
                if len(p) > 1: cl = float(p[-1])
            elif s.startswith("CMytot,"):
                p = s.split(",")
                if len(p) > 1: cm_ = float(p[-1])
    return cl, cm_


def execute_simulation(alpha, delta_e, parm_id, csv_path):
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
    vsp.WriteResultsCSVFile(rid, csv_path)
    vsp.DeleteAllResults()
    return _read_cl_cm(csv_path)


def run_trim(parm_id, CL_objective, csv_path, a0=2.0, de0=0.0):
    """a0/de0: warm-start guess (pass the previous AR point's converged trim
    to cut the Newton iteration down to a handful of steps -- trim varies
    smoothly with AR, so this lands very close to the answer immediately,
    unlike always restarting from a=2 deg, de=0 deg)."""
    cl0, cm0 = execute_simulation(0.0, 0.0, parm_id, csv_path)
    cl2, _   = execute_simulation(2.0, 0.0, parm_id, csv_path)
    _,  cmd5 = execute_simulation(0.0, 5.0, parm_id, csv_path)
    CL_alpha = (cl2 - cl0) / 2.0
    Cm_delta = (cmd5 - cm0) / 5.0
    a, de, trimmed = a0, de0, False
    for _ in range(max_iter):
        CL, CM = execute_simulation(a, de, parm_id, csv_path)
        if CL is None:
            break
        eCL = CL_objective - CL
        eCm = -CM
        if abs(eCL) < tolerance and abs(eCm) < tolerance:
            trimmed = True
            break
        a  += eCL / CL_alpha
        de += eCm / Cm_delta
        a  = max(-4.0, min(a, 10.0))
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


def run_stability(a_trim, de_trim, parm_id, flt_path):
    an = "VSPAEROSweep"
    vsp.SetParmVal(parm_id, de_trim)
    vsp.SetDoubleAnalysisInput(an, "AlphaStart", [a_trim], 0)
    vsp.SetDoubleAnalysisInput(an, "AlphaEnd", [a_trim], 0)
    vsp.SetIntAnalysisInput(an, "AlphaNpts", [1], 0)
    vsp.SetIntAnalysisInput(an, "UnsteadyType", [1], 0)
    vsp.Update()
    vsp.SetAnalysisInputDefaults("VSPAEROComputeGeometry")
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "ThinGeomSet", [0], 0)
    vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", [-1], 0)
    vsp.ExecAnalysis("VSPAEROComputeGeometry")
    vsp.ExecAnalysis(an)
    return parse_flt_file(flt_path, a_trim)


def build_state_matrices(flt, mp, a_trim, S, bref, mac):
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

    Cxu = -(CDu + 2*CD0); Cxa = CL0 - CDa
    Czu = -(CLu + 2*CL0); Cza = -(CLa + CD0); Czq = -CLq

    Xu = 0.5*rho*uo*S*Cxu;  Xw = 0.5*rho*uo*S*Cxa
    Zu = 0.5*rho*uo*S*Czu;  Zw = 0.5*rho*uo*S*Cza;  Zq = 0.25*rho*uo*c*S*Czq
    Mu = 0.5*rho*uo*c*S*Cmu
    Mw = 0.5*rho*uo*c*S*Cma
    Mq = 0.25*rho*uo*c*c*S*Cmq

    A_long = np.array([
        [Xu/m, Xw/m, 0.0,          -g*np.cos(th)],
        [Zu/m, Zw/m, (Zq+m*uo)/m,  -g*np.sin(th)],
        [Mu/Iy, Mw/Iy, Mq/Iy,       0.0],
        [0.0,  0.0,   1.0,          0.0]])

    Yv = 0.5*rho*uo*S*CYb;  Yp = 0.25*rho*uo*S*b*CYp;  Yr = 0.25*rho*uo*S*b*CYr
    Lv = 0.5*rho*uo*S*b*Clb; Lp = 0.25*rho*uo*S*b*b*Clp; Lr = 0.25*rho*uo*S*b*b*Clr
    Nv = 0.5*rho*uo*S*b*Cnb; Np = 0.25*rho*uo*S*b*b*Cnp; Nr = 0.25*rho*uo*S*b*b*Cnr
    Ixp = (Ix*Iz - Ixz**2)/Iz; Izp = (Ix*Iz - Ixz**2)/Ix; Izxp = Ixz/(Ix*Iz - Ixz**2)

    A_lat = np.array([
        [Yv/m,            Yp/m,            (Yr/m)-uo,       g*np.cos(th)],
        [Lv/Ixp+Izxp*Nv, Lp/Ixp+Izxp*Np,  Lr/Ixp+Izxp*Nr, 0.0],
        [Izxp*Lv+Nv/Izp, Izxp*Lp+Np/Izp,  Izxp*Lr+Nr/Izp, 0.0],
        [0.0,            1.0,             np.tan(th),      0.0]])

    SM = -Cma/CLa if CLa != 0 else np.nan
    return A_long, A_lat, SM, CLa, Cma


def _osc(ev):
    out = []
    for e in ev:
        if e.imag > 1e-6:
            wn = abs(e)
            out.append((wn, -e.real/wn, e))
    out.sort(key=lambda t: -t[0])
    return out


def extract_long(ev):
    osc = _osc(ev)
    if len(osc) >= 2:
        sp, ph = osc[0], osc[-1]
    elif len(osc) == 1:
        sp, ph = (np.nan, np.nan, None), osc[0]
    else:
        sp = ph = (np.nan, np.nan, None)
    return sp[0], sp[1], ph[0], ph[1], float(np.max(ev.real))


def extract_lat(ev):
    osc = _osc(ev)
    dr = osc[0] if osc else (np.nan, np.nan, None)
    reals = sorted([e.real for e in ev if abs(e.imag) <= 1e-6])
    roll   = reals[0]  if len(reals) >= 1 else np.nan
    spiral = reals[-1] if len(reals) >= 2 else np.nan
    return dr[0], dr[1], roll, spiral


def _mac(a, b):
    """Modal Assurance Criterion between two eigenvectors: 1 = same mode
    shape, 0 = unrelated shapes. Scale/phase-independent."""
    num = abs(np.vdot(a, b)) ** 2
    den = np.vdot(a, a).real * np.vdot(b, b).real
    return num / den


def _best_pair_score(vecs_a, vecs_b):
    m11 = _mac(vecs_a[0], vecs_b[0]); m22 = _mac(vecs_a[1], vecs_b[1])
    m12 = _mac(vecs_a[0], vecs_b[1]); m21 = _mac(vecs_a[1], vecs_b[0])
    return max(m11 + m22, m12 + m21)


def _group_blocks(ev):
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


def _track_sequence(seq):
    """Original one-directional MAC-continuity tracking (see vspaerotrim_bucle.py):
    seeds SP/phugoid from seq[0] by eigenvalue magnitude, then tracks forward
    through the rest of seq by matching eigenvector shapes (MAC), not just
    eigenvalue position -- survives the crossover near the neutral point.
    Requires seq[0] to be a well-behaved point (stable, clean complex SP +
    phugoid pair) since everything else is matched against it."""
    n = len(seq)
    sp_wn = np.full(n, np.nan); sp_z = np.full(n, np.nan)
    ph_wn = np.full(n, np.nan); ph_z = np.full(n, np.nan)
    sp_ev_list = [None] * n; ph_ev_list = [None] * n

    ev0, evec0 = seq[0]["ev_long"], seq[0]["evec_long"]
    order0 = list(np.argsort(-np.abs(ev0)))
    prev_sp_idx, prev_ph_idx = order0[:2], order0[2:]
    prev_sp_vecs = [evec0[:, k] for k in prev_sp_idx]
    prev_ph_vecs = [evec0[:, k] for k in prev_ph_idx]

    for i, r in enumerate(seq):
        ev, evec = r["ev_long"], r["evec_long"]
        blocks = _group_blocks(ev)
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


def track_longitudinal_modes(results, seed_index=0):
    """Same MAC-continuity tracking as the CG sweep, but anchored at
    `seed_index` instead of always at results[0].

    Why: the original CG sweep is a monotonic sweep that *starts* at a
    stable, well-behaved point, so seeding from index 0 is safe. The AR
    sweep isn't -- with x_cg held fixed, the low-AR end can be genuinely
    unstable (short period already split into two real roots there), so
    seeding from results[0] there poisons the whole tracked sequence (every
    later point gets matched against garbage reference vectors). Instead,
    seed from the known-good baseline point and track outward in both
    directions from it."""
    n = len(results)
    sp_wn = np.full(n, np.nan); sp_z = np.full(n, np.nan)
    ph_wn = np.full(n, np.nan); ph_z = np.full(n, np.nan)
    sp_ev_list = [None] * n; ph_ev_list = [None] * n

    # forward: seed_index, seed_index+1, ..., n-1
    fwd = results[seed_index:]
    f_sp_wn, f_sp_z, f_ph_wn, f_ph_z, f_sp_ev, f_ph_ev = _track_sequence(fwd)
    for j, i in enumerate(range(seed_index, n)):
        sp_wn[i], sp_z[i] = f_sp_wn[j], f_sp_z[j]
        ph_wn[i], ph_z[i] = f_ph_wn[j], f_ph_z[j]
        sp_ev_list[i], ph_ev_list[i] = f_sp_ev[j], f_ph_ev[j]

    # backward: seed_index, seed_index-1, ..., 0 (seed_index already set above)
    if seed_index > 0:
        bwd = results[seed_index::-1]
        b_sp_wn, b_sp_z, b_ph_wn, b_ph_z, b_sp_ev, b_ph_ev = _track_sequence(bwd)
        for j, i in enumerate(range(seed_index, -1, -1)):
            if i == seed_index:
                continue
            sp_wn[i], sp_z[i] = b_sp_wn[j], b_sp_z[j]
            ph_wn[i], ph_z[i] = b_ph_wn[j], b_ph_z[j]
            sp_ev_list[i], ph_ev_list[i] = b_sp_ev[j], b_ph_ev[j]

    return sp_wn, sp_z, ph_wn, ph_z, sp_ev_list, ph_ev_list


# ============================== ASPECT RATIO SWEEP ==============================
AR_mw_values = [3.2, 4.8, 6.4, 8.0, 9.6]   # OpenVSP panel AR -> standard AR (b^2/S) = 2x this
results = []
a0, de0 = 2.0, 0.0   # warm-start guess, updated after each point converges

for AR_mw in AR_mw_values:
    AR_std = 2.0 * AR_mw
    print("\n" + "#"*70)
    print(f"#  AR_mw = {AR_mw:.2f}  (standard AR = b^2/S = {AR_std:.2f})")
    print("#"*70)

    wing_id, htail_id, vtp_id = build_geometry(AR_mw)

    S = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalArea", "WingGeom"))
    bref = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalSpan", "WingGeom"))
    mac = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalChord", "WingGeom"))
    q_dyn = 0.5 * rho * v**2
    CL_objective = W / (q_dyn * S)
    re_cref = rho * v * mac / mu
    print(f"   S={S:.4f} m^2  bref={bref:.4f} m  mac={mac:.4f} m  CL_obj={CL_objective:.4f}")

    mp = build_mass_model(X_CG_FIXED, AR_std)
    print(f"   m={mp['m']:.3f} kg  CG_x={mp['cg'][0]:.4f}  "
          f"Ixx={mp['I1'][0]:.3f} Iyy={mp['I1'][1]:.3f} Izz={mp['I1'][2]:.3f} Ixz={mp['I2'][1]:.3f}")

    vsp.ClearVSPModel(); vsp.ReadVSPFile(BASE_VSP); vsp.Update()
    parm_id = setup_vspaero_inputs(X_CG_FIXED, S, bref, mac, mach, re_cref)

    csv_path = os.path.join(HERE, "Aircraft_Stability.csv")
    flt_path = os.path.join(HERE, "AC.flt")

    a_trim, de_trim, CLa_fit, Cmd_fit, trimmed = run_trim(parm_id, CL_objective, csv_path, a0=a0, de0=de0)
    if not trimmed:
        print(f"   [!] trim not converged at AR={AR_std:.2f} -> skipped")
        continue
    print(f"   TRIM: alpha={a_trim:.4f} deg | elevator={de_trim:.4f} deg")
    a0, de0 = a_trim, de_trim   # warm-start the next AR point from here

    flt = run_stability(a_trim, de_trim, parm_id, flt_path)
    A_long, A_lat, SM, CLa, Cma = build_state_matrices(flt, mp, a_trim, S, bref, mac)
    ev_long, evec_long = np.linalg.eig(A_long)
    ev_lat = np.linalg.eigvals(A_lat)

    dr_wn, dr_z, roll, spiral = extract_lat(ev_lat)

    print(f"   SM={SM*100:6.2f}%MAC")

    results.append(dict(AR_mw=AR_mw, AR_std=AR_std, S=S, bref=bref, mac=mac,
                        alpha=a_trim, dele=de_trim, SM=SM,
                        Ixx=mp["I1"][0], Iyy=mp["I1"][1], Izz=mp["I1"][2], Ixz=mp["I2"][1],
                        ev_long=ev_long, evec_long=evec_long, ev_lat=ev_lat,
                        dr_wn=dr_wn, dr_z=dr_z, roll=roll, spiral=spiral))

# ------ track short period / phugoid by continuity across the AR sweep ------
# Anchor the tracking at the baseline AR (AR_mw=6.4, the original design point,
# stable at x_cg=1.325m) rather than always at results[0] -- see
# track_longitudinal_modes docstring for why: with x_cg held fixed, the
# low-AR end of this particular sweep is unstable, so seeding from there
# would poison the tracking for every other point.
_seed_i = min(range(len(results)), key=lambda k: abs(results[k]["AR_mw"] - 6.4))
print(f"\nmode tracking anchored at AR_mw={results[_seed_i]['AR_mw']:.2f} "
      f"(AR_std={results[_seed_i]['AR_std']:.2f}, SM={results[_seed_i]['SM']*100:.2f}%MAC)")
sp_wn_trk, sp_z_trk, ph_wn_trk, ph_z_trk, sp_ev_list, ph_ev_list = track_longitudinal_modes(results, seed_index=_seed_i)

# Per-point check only (no cascading "once collapsed, stay collapsed"): that
# cascade rule made sense for the CG sweep (a monotonic sweep that starts
# stable and only ever crosses into instability once, aft of the neutral
# point). Here the sweep is anchored at a stable baseline and radiates both
# ways, so an unstable point at one end must not blank out perfectly good
# points further along -- each point is judged on its own eigenvalues only.
for i, r in enumerate(results):
    r["collapsed"] = len(_osc(r["ev_long"])) == 0

for i, r in enumerate(results):
    r["sp_wn"], r["sp_z"] = sp_wn_trk[i], sp_z_trk[i]
    r["ph_wn"], r["ph_z"] = ph_wn_trk[i], ph_z_trk[i]
    r["sp_ev"], r["ph_ev"] = sp_ev_list[i], ph_ev_list[i]

# ============================== SAVE TABLE ==============================
out_csv = os.path.join(HERE, "ar_sweep_aero_results.csv")
with open(out_csv, "w") as f:
    f.write("AR_mw,AR_standard,S,bref,mac,alpha_trim,elevator_trim,SM_%MAC,Iyy,"
            "SP_freq_Hz,SP_zeta,Phug_freq_Hz,Phug_zeta,DR_freq_Hz,DR_zeta,roll_eig,spiral_eig\n")
    for r in results:
        f.write(f"{r['AR_mw']:.3f},{r['AR_std']:.3f},{r['S']:.4f},{r['bref']:.4f},{r['mac']:.4f},"
                f"{r['alpha']:.4f},{r['dele']:.4f},{r['SM']*100:.3f},{r['Iyy']:.4f},"
                f"{r['sp_wn']/(2*np.pi):.4f},{r['sp_z']:.4f},"
                f"{r['ph_wn']/(2*np.pi):.4f},{r['ph_z']:.4f},"
                f"{r['dr_wn']/(2*np.pi):.4f},{r['dr_z']:.4f},{r['roll']:.4f},{r['spiral']:.4f}\n")
print(f"\nresults table saved: {out_csv}")

# ============================== PLOTS ==============================
arr = {k: np.array([r[k] for r in results]) for k in
       ("AR_std", "alpha", "dele", "SM", "Iyy", "sp_wn", "sp_z", "ph_wn", "ph_z", "bref", "mac")}

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(arr["AR_std"], arr["sp_wn"]/(2*np.pi), "o-", color="tab:blue")
ax.set_xlabel("Aspect ratio (b$^2$/S)"); ax.set_ylabel("Short-period frequency [Hz]")
ax.set_title("Short-period frequency vs. aspect ratio (S fixed at 0.8 m$^2$)")
ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "ar_sweep_short_period.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(arr["AR_std"], arr["bref"], "o-", color="tab:green", label="span")
ax.plot(arr["AR_std"], arr["mac"], "s--", color="tab:orange", label="MAC")
ax.set_xlabel("Aspect ratio (b$^2$/S)"); ax.set_ylabel("Length [m]")
ax.set_title("Span and mean aerodynamic chord vs. aspect ratio")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "ar_sweep_geometry.png"), dpi=150); plt.close(fig)

print("\n" + "="*60)
print("AR SWEEP (aero) COMPLETE")
print(f"  cases run: {len(results)}")
print(f"  data : {out_csv}")
print("="*60)
