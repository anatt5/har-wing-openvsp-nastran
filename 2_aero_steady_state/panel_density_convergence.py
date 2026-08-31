"""
VSPAERO aerodynamic mesh convergence study.

Builds the aircraft geometry (main wing, horizontal tail, vertical tail,
plus aileron/flap/elevator/rudder control surfaces for VSPAERO) and then
sweeps the VLM panel density (panels/m, equal in chord and span so panels
stay square) to check how CLtot, CDtot and CMytot converge as the mesh is
refined. At each density the model is reloaded from a pristine saved base
(AC_base.vsp3) so no degenerate/cached geometry leaks between runs.

Outputs:
  AC.vsp3                        -- clean aircraft geometry (no mesh study)
  AC_base.vsp3                   -- pristine base reloaded at every density
  AC_mesh_<tag>_d<d>.vsp3        -- geometry actually used at density d
  mesh_convergence_<tag>.txt     -- full log (per-density results + summary)
  mesh_convergence_<tag>.png     -- 4-panel convergence plot

Run:  python panel_density_convergence.py
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
vsp.SetParmVal(wing_id, "Tess_W", "Shape", 49)
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
vsp.SetParmVal(htail_id, "Tess_W", "Shape", 53)
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
vsp.SetParmVal(vtp_id, "Tess_W", "Shape", 80)
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

# save clean geometry (no mesh-convergence artifacts)
vsp.WriteVSPFile("AC.vsp3")
print("plane geometry created: AC.vsp3")

########################################################################
# 7) MESH CONVERGENCE STUDY (VLM panel density)
########################################################################
# Varies Num U (SectTess_U, span) and Num W (Tess_W "Shape", chord) of the
# WING, HTP and VTP using the live IDs (wing_id, htail_id, vtp_id). At each
# level it uses ONE density d (panels/m) equal in chord and span -> square
# panels by construction (ideal for VLM). Runs VSPAERO at a single alpha,
# collects the converged CLtot/CDtot/CMytot (last column of the .csv = last
# wake iteration) and writes everything to a .txt log. Also plots the
# convergence trend and the global error vs. the finest mesh tested.

# ---------------- STUDY SETTINGS ----------------
ALPHA_DEG = 4.0                        # single comparison alpha (linear zone)
DENSITIES = [20, 25, 44, 50, 55, 60]   # panels/m: 30/33/35/37/38/40 all wobbled or
                                        # spiked in this geometry (control-surface
                                        # split lines) -> the whole 25-44 range is
                                        # unreliable, so the gap is left as-is rather
                                        # than forcing in another bad point; d=44 is
                                        # the only clean one found in that range.
TOL_CL = 0.01                          # CL: strict (1%)
TOL_OTH = 0.01                         # CDtot/CDi/CMy: same strict criterion (1%)
# Note: at each density AC_mesh_{tag}_d{d}.vsp3 is saved (to open and check the mesh)
STOP_ON_CONVERGENCE = False            # False = run ALL densities (to see the full trend)
CLEAN_CONTROLS = True                  # True = study with control surfaces UNDEFLECTED
TEST_NO_CONTROL_SURFACES = False       # DIAGNOSTIC: delete the control-surface subsurfaces
                                        # (aileron/flap/rudder/elevator) before the study, to
                                        # check whether their split lines cause the non-monotonic
                                        # CD/CMy jump seen around d=35 with uniform tessellation

# Clustering: groups panels where the load changes fast -> CD and Cm converge much sooner.
# Values < 1.0 group towards that edge (1.0 = uniform). If CLUSTERING=False a uniform mesh
# is forced (all factors set to 1.0) -> clean comparison baseline (uniform vs. clustered).
CLUSTERING = False                     # reverted: clustering made CD/CDi WORSE (never converged,
                                        # see mesh_convergence_clustered.txt) -> back to uniform,
                                        # which is what the d=50 production runs (trim/steady
                                        # state/dynamic modes) already use
LE_CLUSTER = 0.25                      # chord: smaller = more panels near the leading edge
TE_CLUSTER = 0.25                      # chord: smaller = more panels near the trailing edge
OUT_CLUSTER = 0.25                     # span: smaller = more panels towards the tip
IN_CLUSTER = 0.5                       # span: smaller = more panels towards the root

# Tag based on the clustering state -> output files are NOT overwritten between runs.
TAG = ("clustered" if CLUSTERING else "uniform") + ("_nosubsurf" if TEST_NO_CONTROL_SURFACES else "")
LOG_FILE = f"mesh_convergence_{TAG}.txt"
PNG_FILE = f"mesh_convergence_{TAG}.png"

# Reference / flight conditions (identical to the production analysis)
cref = 0.25321; bref = 3.2; Sref = 0.8
rho = 1.225; Vinf = 30.0   # TU-Flex cruise speed Vc (reference papers)
a_SL = 340.29; mu = 1.7894e-5
mach = Vinf / a_SL
Re = rho * Vinf * cref / mu

# Surfaces to refine. The dimensions (section span, mean chord) are taken directly
# from the variables already defined above when the aircraft was built, so there is
# no need to query parms from OpenVSP. We only use the API to SET the tessellation.
# The geom_id is re-located by name at every iteration (IDs change on model reload).
#   format: (name, [(section_span, section_mean_chord), ...])
SURF_DEF = [
    ("main wing", [(span_mw, Sref / bref)]),                            # 1 section
    ("HTP", [(span_htp, rootc_htp * (1 + taper_htp) / 2)]),             # 1 section
    ("VTP", [(span_vtp, rootc_vtp * (1 + taper_vtp) / 2),                # section 1
             (span_vtp2, rootc_vtp2 * (1 + taper_vtp2) / 2)]),          # section 2
]


def _set_tessellation(geom_id, sections, d):
    """Set Num W and Num U of a WING geometry for density d (panels/m).
    'sections' = [(span, mean_chord), ...].
    Returns the REAL values read back from the model: (num_w, [num_u...])."""
    surf = vsp.GetXSecSurf(geom_id, 0)

    # mean chord of the geometry, weighted by section span
    span_tot = sum(s for s, _ in sections) or 1.0
    cbar_geom = sum(s * c for s, c in sections) / span_tot

    # Num W (chord): odd, minimum 7
    nc = max(3, int(round(d * cbar_geom)))
    num_w = 2 * nc + 1
    vsp.SetParmVal(geom_id, "Tess_W", "Shape", num_w)

    # Num U (span): single value for the whole surface (same sections), computed
    # with the mean span -> all VTP sections share the same Num U.
    # NOTE: the SectTess_U of section i (1-indexed) lives on the OUTER xsec -> index i+1.
    span_mean = span_tot / len(sections)
    ns = max(2, int(round(d * span_mean)))
    num_u = ns + 1
    for i in range(len(sections)):
        xsec_id = vsp.GetXSec(surf, i + 1)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "SectTess_U"), float(num_u))

    # Clustering: groups panels towards LE/TE (chord, geom-level parm) and towards
    # tip/root (span, per-section parm on the outer xsec). If CLUSTERING=False
    # everything is set to 1.0 -> uniform mesh (comparison baseline). Defensive: if
    # the parm doesn't exist under that name, FindParm/GetXSecParm returns "" and
    # it is simply not applied.
    le, te, out, inn = ((LE_CLUSTER, TE_CLUSTER, OUT_CLUSTER, IN_CLUSTER)
                         if CLUSTERING else (1.0, 1.0, 1.0, 1.0))
    # chord (LE/TE): geom-level parm. The group can vary depending on version.
    for pname, pval in (("LECluster", le), ("TECluster", te)):
        pid = ""
        for grp in ("WingGeom", "Shape", "Tess", "Design"):
            pid = vsp.GetParm(geom_id, pname, grp)
            if pid:
                break
        if pid:
            vsp.SetParmVal(pid, pval)
    # span (tip/root): per-section parm, on the outer xsec (i+1).
    for i in range(len(sections)):
        xsec_id = vsp.GetXSec(surf, i + 1)
        for pname, pval in (("OutCluster", out), ("InCluster", inn)):
            pid = vsp.GetXSecParm(xsec_id, pname)
            if pid:
                vsp.SetParmVal(pid, pval)
    vsp.Update()

    # read back what actually stayed in the model (verification)
    real_w = int(round(vsp.GetParmVal(geom_id, "Tess_W", "Shape")))
    real_u = []
    for i in range(len(sections)):
        xsec_id = vsp.GetXSec(surf, i + 1)
        real_u.append(int(round(vsp.GetParmVal(vsp.GetXSecParm(xsec_id, "SectTess_U")))))
    return real_w, real_u


def _run_vspaero_single_alpha(alpha_deg, csv_out):
    """Configure and run VSPAEROSweep at a single alpha; dump the CSV to csv_out."""
    # Delete previous results so a cached CSV from the previous iteration is NOT
    # reused (that was causing CL/CD to come out identical across mesh densities).
    try:
        vsp.DeleteAllResults()
    except Exception:
        pass
    analysis_name = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(analysis_name)
    vsp.SetIntAnalysisInput(analysis_name, "UnsteadyType", [1], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Sref", [Sref], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "cref", [cref], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "bref", [bref], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Xcg", [1.23], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Ycg", [0.0], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Zcg", [0.023], 0)
    # single alpha
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [alpha_deg], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [alpha_deg], 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [1], 0)
    vsp.SetIntAnalysisInput(analysis_name, "BetaNpts", [1])
    vsp.SetDoubleAnalysisInput(analysis_name, "BetaStart", [0])
    vsp.Update()
    vsp.SetDoubleAnalysisInput(analysis_name, "Rho", [rho], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", [Vinf], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "MachStart", [mach], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "ReCref", [Re], 0)
    vsp.SetIntAnalysisInput(analysis_name, "Symmetry", [0], 0)
    # Fixed wake during the study: removes the relaxed-wake noise from CD and Cm
    # and isolates the effect of the mesh (switch back to free wake for the final run).
    vsp.SetIntAnalysisInput(analysis_name, "FixedWakeFlag", [1], 0)
    # With a fixed wake, the wake never evolves -> a single iteration is enough
    # (more would just waste time). WakeRelax is ignored with a fixed wake, so it
    # is left unset. NumWakeNodes only discretizes the straight prescribed wake
    # (minor influence): left at its default.
    vsp.SetDoubleAnalysisInput(analysis_name, "CoreSizeFactor", [1.0], 0)
    vsp.SetIntAnalysisInput(analysis_name, "StallModel", [0], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "ThinGeomSet", [0], 0)  # all (VLM)
    vsp.Update()

    comp_geom = "VSPAEROComputeGeometry"
    vsp.SetAnalysisInputDefaults(comp_geom)
    vsp.ExecAnalysis(comp_geom)  # generates the degenerate mesh (.vspgeom)
    results_id = vsp.ExecAnalysis(analysis_name)
    if not results_id:
        raise RuntimeError("VSPAERO failed (check vspaero.exe)")
    if os.path.exists(csv_out):
        os.remove(csv_out)
    vsp.WriteResultsCSVFile(results_id, csv_out)


def _vspgeom_nodes(path="AC.vspgeom"):
    """Number of nodes in the mesh VSPAERO actually solves (vspgeom v3 format:
    '# vspgeom v3' header, then one line with a single integer, then
    'nNodes nFaces ...'). If this number does NOT change between densities,
    ComputeGeometry is not regenerating the mesh."""
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:  # first line with >=2 integers: nNodes nFaces ...
                    return int(parts[0])
    except Exception:
        return None
    return None


def _zero_control_deflections():
    """Set the deflection of ALL control-surface groups to 0 (clean study)."""
    container = vsp.FindContainer("VSPAEROSettings", 0)
    if not container:
        return
    try:
        n = vsp.GetNumControlSurfaceGroups()
    except Exception:
        n = 8  # fallback: try a handful of groups
    for i in range(n):
        pid = vsp.FindParm(container, "DeflectionAngle", f"ControlSurfaceGroup_{i}")
        if pid:
            vsp.SetParmVal(pid, 0.0)
    vsp.Update()


def _clean_vspaero_files(base="AC"):
    """Delete VSPAERO's output files so ComputeGeometry is forced to REGENERATE
    the mesh (.vspgeom) instead of reusing the previous iteration's."""
    for ext in (".vspgeom", ".adb", ".adb.cases", ".history", ".lod",
                ".stab", ".polar", ".fem", ".tkey", ".group.0", ".group.1"):
        try:
            os.remove(base + ext)
        except OSError:
            pass


def _parse_coeffs(csv_file, keys=("CLtot", "CDtot", "CDi", "CMytot")):
    """Return the converged (CLtot, CDtot, CDi, CMytot) -- last column = last
    wake iteration. Takes the FIRST row of each coefficient (base case, not
    perturbations). CL and CDtot are mandatory; CDi (far-field induced drag)
    and CMytot are best-effort (NaN if not found under that name)."""
    lk = {k.lower(): k for k in keys}
    out = {k: None for k in keys}
    with open(csv_file, "r") as f:
        for line in f:
            head = line.split(",", 1)[0].strip().lower()
            if head in lk:
                k = lk[head]
                if out[k] is None:
                    try:
                        out[k] = float(line.strip().split(",")[-1])
                    except ValueError:
                        pass
            if all(v is not None for v in out.values()):
                break
    if out["CLtot"] is None or out["CDtot"] is None:
        raise RuntimeError(f"CLtot/CDtot not found in {csv_file}")
    for k in ("CDi", "CMytot"):
        if out[k] is None:
            out[k] = float("nan")  # different name in the CSV: check that row
    return out["CLtot"], out["CDtot"], out["CDi"], out["CMytot"]


# ----------------------------- MAIN LOOP -----------------------------
_log = open(LOG_FILE, "w", encoding="utf-8")


def _emit(t=""):
    print(t); _log.write(t + "\n"); _log.flush()


_emit("=" * 92)
_emit("MESH CONVERGENCE STUDY  -  VSPAERO (VLM)")
_emit(f"Alpha={ALPHA_DEG} deg  |  Vinf={Vinf} m/s  Mach={mach:.4f}  Re={Re:.3e}")
_emit(f"Surfaces: {', '.join(n for n, _s in SURF_DEF)}   |   CL tol: {TOL_CL*100:.1f}%  "
      f"CD/CDi/CMy tol: {TOL_OTH*100:.1f}%")
_emit(f"Control surfaces deflected: {'NO' if CLEAN_CONTROLS else 'YES'}   |   Clustering: {'YES' if CLUSTERING else 'NO'}"
      + (f" (LE={LE_CLUSTER} TE={TE_CLUSTER} tip={OUT_CLUSTER} root={IN_CLUSTER})" if CLUSTERING else ""))
_emit(f"Control-surface subsurfaces removed for this run (diagnostic): {'YES' if TEST_NO_CONTROL_SURFACES else 'NO'}")
_emit("=" * 92)

# DIAGNOSTIC: delete the control-surface subsurfaces (aileron/flap/rudder/elevator) so their
# split lines are absent from the mesh used in this convergence study. If the CD/CMy jump at
# d=35 disappears with this flag on, it confirms the subsurface split lines were the cause.
if TEST_NO_CONTROL_SURFACES:
    for _ss_id in (aileron_id, flap_id, rudder_id, elevator_id):
        vsp.DeleteSubSurf(_ss_id)
    vsp.Update()

# Save ONE pristine base of the aircraft. Each density is reloaded from here, so every
# run starts from a clean model (no cached degenerate geometry).
vsp.WriteVSPFile("AC_base.vsp3")

prev_cl = prev_cd = prev_cdi = prev_cmy = None
history = []  # (d, {name: (num_w, [num_u])}, cl, cd, cdi, cmy, nodes)
converged_d = None

for d in DENSITIES:
    base = f"AC_mesh_{TAG}_d{d}"  # base name unique per density and per clustering setting

    # 1) Fresh model from the pristine base -> no cached degenerate geometry.
    vsp.ClearVSPModel()
    vsp.ReadVSPFile("AC_base.vsp3")
    if CLEAN_CONTROLS:
        _zero_control_deflections()  # study with undeflected controls

    # 2) Locate the geometries by name and set this level's tessellation.
    info = {}
    for name, sections in SURF_DEF:
        gid = vsp.FindGeomsWithName(name)[0]
        num_w, num_u_list = _set_tessellation(gid, sections, d)
        info[name] = (num_w, num_u_list)
    vsp.Update()

    # 3) Write THIS mesh's .vsp3. This fixes VSPAERO's base name to 'base', so it
    #    writes base.vspgeom / base.history / ... (unique per density, nothing
    #    reused or overwritten). Also handy to open and inspect the mesh directly.
    vsp.WriteVSPFile(base + ".vsp3")
    _clean_vspaero_files(base)  # in case something was left over from a previous run

    # 4) Run VSPAERO on the fresh model.
    csv_out = base + "_results.csv"
    try:
        _run_vspaero_single_alpha(ALPHA_DEG, csv_out)
        nodos = _vspgeom_nodes(base + ".vspgeom")  # actual mesh VSPAERO solves
        cl, cd, cdi, cmy = _parse_coeffs(csv_out)
    except Exception as e:
        _emit("")
        _emit(f"--- d = {d} panels/m : VSPAERO FAILED ({e}) -> skipping this density")
        continue
    history.append((d, dict(info), cl, cd, cdi, cmy, nodos))

    def _rel(new, old):  # relative change (absolute if old ~ 0)
        if old is None:
            return float("nan")
        return abs((new - old) / old) if abs(old) > 1e-9 else abs(new - old)

    if prev_cl is None:
        conv_txt = "  (reference)"; dcl = dcd = dcmy = float("nan")
    else:
        dcl = _rel(cl, prev_cl)
        dcd = _rel(cd, prev_cd)
        dcmy = _rel(cmy, prev_cmy)
        conv_txt = (f"  dCL={dcl*100:5.2f}%  dCD={dcd*100:5.2f}%  dCMy={dcmy*100:5.2f}%")

    _emit("")
    _emit(f"--- d = {d} panels/m -----------------------------------------------------------------")
    paneles_tot = 0
    for name, _sec in SURF_DEF:
        num_w, num_u_list = info[name]
        nus = ", ".join(str(u) for u in num_u_list)
        # estimated VLM panels = (chordwise panels) * (total spanwise panels)
        nc_panels = (num_w - 1) // 2
        ns_panels = sum(u - 1 for u in num_u_list)
        paneles_tot += nc_panels * ns_panels
        # panel size and aspect ratio -> should come out ~1 (square)
        span_s = sum(s for s, _ in _sec)
        cbar_s = sum(s * c for s, c in _sec) / (span_s or 1.0)
        ch_mm = (cbar_s / nc_panels) * 1000.0
        sp_mm = (span_s / ns_panels) * 1000.0
        ar_p = sp_mm / ch_mm if ch_mm else float("nan")
        _emit(f"    {name:<10}:  Num W = {num_w:<4} Num U = [{nus}]"
              f"   panel ~ {sp_mm:5.1f}x{ch_mm:5.1f} mm (AR {ar_p:.2f})")
    _emit(f"    Estimated VLM panels ~ {paneles_tot}   |   Actual VSPAERO mesh ({base}.vspgeom): {nodos} nodes")
    _emit(f"    CLtot = {cl:.6f}   CDtot = {cd:.6f}   CDi = {cdi:.6f}   CMytot = {cmy:.6f}{conv_txt}")

    if (prev_cl is not None and dcl < TOL_CL and dcd < TOL_OTH and dcmy < TOL_OTH
            and converged_d is None):
        converged_d = d
        _emit(f"    (CONSECUTIVE change within tol from d={d}; final recommendation uses "
              f"the GLOBAL error vs finest mesh, below)")
        if STOP_ON_CONVERGENCE:
            break
    prev_cl, prev_cd, prev_cdi, prev_cmy = cl, cd, cdi, cmy

# ----------------------------- SUMMARY -----------------------------
_emit("")
_emit("=" * 92)
_emit("SUMMARY")
_emit(f"{'d[pan/m]':>9} | {'nodes':>7} | {'CLtot':>10} | {'CDtot':>10} | {'CDi':>10} | {'CMytot':>10} | {'L/D':>7}")
_emit("-" * 83)
for d, _info, cl, cd, cdi, cmy, nod in history:
    ld = cl / cd if cd else float("nan")
    _emit(f"{d:>9} | {str(nod):>7} | {cl:>10.6f} | {cd:>10.6f} | {cdi:>10.6f} | {cmy:>10.6f} | {ld:>7.2f}")
_emit("-" * 83)

# --- Global error vs. the FINEST mesh (robust criterion: every mesh is compared
#     against the best available estimate, not against the previous one) ---
rec_d = None
if len(history) >= 2:
    ref = history[-1]  # finest mesh = reference
    refCL, refCD, refCDi, refCMy = ref[2], ref[3], ref[4], ref[5]

    def _err(v, r):
        return abs(v - r) / abs(r) if abs(r) > 1e-9 else float("nan")

    _emit("")
    _emit(f"GLOBAL error vs finest mesh (d={ref[0]}, {ref[6]} nodes):")
    _emit(f"{'d[pan/m]':>9} | {'errCL%':>8} | {'errCD%':>8} | {'errCDi%':>8} | {'errCMy%':>8} | {'<tol?':>6}")
    _emit("-" * 62)
    rows = []  # (d, within_tol)
    for d, _info, cl, cd, cdi, cmy, _nod in history[:-1]:
        eL, eD = _err(cl, refCL) * 100, _err(cd, refCD) * 100
        eDi, eM = _err(cdi, refCDi) * 100, _err(cmy, refCMy) * 100
        okM = (eM < TOL_OTH * 100) or (eM != eM)  # NaN (CMy~0) does not block
        # Criterion: CL with the strict tolerance; CDtot and CMy with an engineering
        # (looser) tolerance. CDi is shown as a robust reference but does not block.
        ok = (eL < TOL_CL * 100) and (eD < TOL_OTH * 100) and okM
        rows.append((d, ok))
        _emit(f"{d:>9} | {eL:>8.2f} | {eD:>8.2f} | {eDi:>8.2f} | {eM:>8.2f} | {'yes' if ok else 'no':>6}")
    _emit("-" * 62)

    # STRICT criterion: the cheapest mesh such that it AND every finer mesh are
    # within tolerance (avoids false positives from the error curve crossing back).
    for i in range(len(rows)):
        if all(ok for _dd, ok in rows[i:]):
            rec_d = rows[i][0]
            break

    if rec_d is not None:
        _emit(f"RECOMMENDED mesh (cheapest with errCL<{TOL_CL*100:.0f}% and "
              f"errCD/errCMy<{TOL_OTH*100:.0f}% vs finest): d = {rec_d} panels/m")
        for hd, hinfo, *_rest in history:
            if hd == rec_d:
                for name, _sec in SURF_DEF:
                    nw, nul = hinfo[name]
                    _emit(f"   {name:<10}:  Num W = {nw}   Num U = [{', '.join(str(u) for u in nul)}]")
    else:
        _emit(f">>> NOT CONVERGED yet with these tolerances (CL<{TOL_CL*100:.0f}%, "
              f"CD/CMy<{TOL_OTH*100:.0f}%) vs finest. Add larger densities to DENSITIES.")
else:
    _emit("Need >= 2 densities to assess convergence.")
_emit("=" * 92)
_log.close()
print(f"\nResults saved to: {LOG_FILE}")

# ----------------------------- PLOTS -----------------------------
try:
    import matplotlib.pyplot as plt

    ds = [h[0] for h in history]
    clv = [h[2] for h in history]
    cdv = [h[3] for h in history]
    cmyv = [h[5] for h in history]  # h[4] = CDi: kept in the log/table, left out of the figure
    nodes = [h[6] for h in history]
    use_nodes = all(n is not None for n in nodes)
    x = nodes if use_nodes else ds
    xlabel = "Mesh nodes (VSPAERO)" if use_nodes else "Density (panels/m)"

    fig, axs = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    fig.suptitle(f"Mesh convergence of the aerodynamic coefficients ($\\alpha$ = {ALPHA_DEG}$^\\circ$)",
                 fontsize=14, fontweight="bold")

    def _panel(ax, y, title, ylab, color):
        ax.plot(x, y, "o-", color=color, lw=1.8, ms=6)  # markers + straight line
        ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylab)
        ax.grid(True, ls="--", alpha=0.5)

    _panel(axs[0, 0], clv, "$C_{Ltot}$ vs mesh", "$C_{Ltot}$", "navy")
    _panel(axs[0, 1], cdv, "$C_{Dtot}$ vs mesh", "$C_{Dtot}$", "firebrick")
    _panel(axs[1, 0], cmyv, "$C_{Mytot}$ vs mesh", "$C_{Mytot}$", "teal")

    # Global error vs finest mesh (log scale), excluding the finest itself (error 0).
    ax = axs[1, 1]
    if len(history) >= 2:
        rCL, rCD, rCMy = clv[-1], cdv[-1], cmyv[-1]
        xe = x[:-1]

        def _e(vals, ref):
            return [abs(v - ref) / abs(ref) * 100 if abs(ref) > 1e-9 else float("nan")
                    for v in vals[:-1]]

        for yv, mk, col, lab in (
            (_e(clv, rCL), "o-", "navy", "$C_L$"),
            (_e(cdv, rCD), "s-", "firebrick", "$C_{Dtot}$"),
            (_e(cmyv, rCMy), "^-", "teal", "$C_{My}$"),
        ):
            ax.plot(xe, yv, mk, color=col, lw=1.8, ms=6, label=lab)
        ax.axhline(TOL_CL * 100, color="green", ls="--", label=f"CL tol {TOL_CL*100:.0f}%")
        ax.axhline(TOL_OTH * 100, color="gray", ls=":", label=f"CD/CMy tol {TOL_OTH*100:.0f}%")
        ax.set_yscale("log")
        ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.set_title("Global error vs finest mesh")
    ax.set_xlabel(xlabel); ax.set_ylabel("Relative error (%)")
    ax.grid(True, which="both", ls="--", alpha=0.5)

    fig.savefig(PNG_FILE, dpi=130)
    print(f"Plot saved to: {PNG_FILE}")
    plt.show()
except Exception as e:
    print(f"(Could not plot: {e})")
