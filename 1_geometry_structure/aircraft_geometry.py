"""
FULL AIRCRAFT GEOMETRY -- wing + horizontal tail + vertical tail + fuselage.

Geometry-only build (no VSPAERO analysis here): loads the "TransportFuse"
template as the fuselage, re-scales it to this aircraft's dimensions, and
builds the main wing, HTP and 2-section VTP on top of it, with their
control surfaces and VSPAERO control-surface groups already defined.
Ends by writing out "AC.vsp3".

Of the three fuselage approaches tried during development (hand-built
elliptical cross-sections, a refined version of the same, and loading the
stock "TransportFuse0.vsp3" template and rescaling it), this script keeps
only the TransportFuse version -- the one actually used going forward.

Needs "TransportFuse0.vsp3" (OpenVSP's stock transport-fuselage example
file) in the same folder as this script.

Run:  python aircraft_geometry.py
"""
import os
import sys

import openvsp as vsp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
AIRFOIL_DIR = os.path.join(REPO_ROOT, "airfoils")
FUSELAGE_TEMPLATE = os.path.join(HERE, "TransportFuse0.vsp3")

sys.path.insert(0, REPO_ROOT)
from utils.airfoil_utils import assign_airfoil_to_component

vsp.VSPRenew()
vsp.VSPCheckSetup()

# ============================================================ FUSELAGE ============================================================
# Load the stock TransportFuse template and rescale it to this aircraft's
# dimensions (found by parameter NAME, not by index, so this survives
# whatever internal cross-section layout the template file uses).
vsp.ReadVSPFile(FUSELAGE_TEMPLATE)

geoms = vsp.FindGeoms()
if geoms:
    fuselage_id = geoms[0]
    print(f"Targeting: {vsp.GetGeomName(fuselage_id)}")

    targets = {
        "Length": 2.7,
        "NoseMult": 1.1,
        "AftMult": 2.95,
        "Diameter": 0.35,
        "NoseCenter": -0.12,
        "AftCenter": 0.33,
        "AftWidth": 0.3,
        "AftHeight": 0.3,
        "Tess_W": 21,
        "Tess_U": 37,
    }

    found_ids = {}
    parm_ids = vsp.GetGeomParmIDs(fuselage_id)
    for p_id in parm_ids:
        p_name = vsp.GetParmName(p_id)
        if p_name in targets:
            found_ids[p_name] = p_id

    for name, value in targets.items():
        if name in found_ids:
            vsp.SetParmVal(found_ids[name], value)
            print(f"Successfully set {name} to {value}")
        else:
            print(f"Warning: Parameter '{name}' not found in this component.")

    vsp.UpdateGeom(fuselage_id)
    print("fuselage finished")

vsp.SetParmVal(fuselage_id, "Density", "Mass_Props", 80)


# ============================================================ MAIN WING ============================================================
wing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(wing_id, "main wing")
# Tess_W=25 is the VSPAERO mesh-convergence result (d=50 panels/m,
# see "AC VSPAERO/Alphasweep at stall 0/mesh convergence/mesh_convergence_uniform.txt"):
# cheapest density with <1% error vs the finest mesh tested.
vsp.SetParmVal(wing_id, "Tess_W", "Shape", 25)

span_mw = 1.6
AR_mw = 6.4
taper_mw = 0.66667
sweep_angle = 28

vsp.SetDriverGroup(wing_id, 1, vsp.SPAN_WSECT_DRIVER, vsp.AR_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
vsp.Update()

vsp.SetParmVal(vsp.GetParm(wing_id, "Span", "XSec_1"), span_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Aspect", "XSec_1"), AR_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Taper", "XSec_1"), taper_mw)
vsp.SetParmVal(vsp.GetParm(wing_id, "Sweep", "XSec_1"), sweep_angle)
vsp.Update()

vsp.SetParmVal(vsp.GetParm(wing_id, "X_Rel_Location", "XForm"), 0.7)
vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Location", "XForm"), 0)
vsp.SetParmVal(vsp.GetParm(wing_id, "Z_Rel_Location", "XForm"), 0)
vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Rotation", "XForm"), 0)
vsp.Update()

assign_airfoil_to_component(wing_id, "clark Y.dat", AIRFOIL_DIR)
xsec_surf_id_mwing = vsp.GetXSecSurf(wing_id, 0)

# NOTE: the mesh-convergence density sweep only ever varies the OUTBOARD
# xsec's SectTess_U (index 1) -- the root xsec (index 0) keeps its own,
# separately-chosen baseline value. Both are carried over as-is from that
# study so this is the exact mesh validated there, not a re-derivation.
xsec_id_1_mwing = vsp.GetXSec(xsec_surf_id_mwing, 0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_1_mwing, "SectTess_U"), 86.0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_1_mwing, "Dihedral"), 6.0)

xsec_id_2_mwing = vsp.GetXSec(xsec_surf_id_mwing, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_2_mwing, "SectTess_U"), 81.0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_2_mwing, "Dihedral"), 6.0)

vsp.SetParmVal(wing_id, "Density", "Mass_Props", 40)
vsp.Update()


# ============================================================ HORIZONTAL TAIL (HTP) ============================================================
htail_id = vsp.AddGeom("WING")
vsp.SetGeomName(htail_id, "HTP")
vsp.SetParmVal(htail_id, "Tess_W", "Shape", 29)  # d=50 mesh-convergence result

span_htp = 0.38
rootc_htp = 0.35
taper_htp = 0.57143

vsp.SetDriverGroup(htail_id, 1, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
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
vsp.Update()

assign_airfoil_to_component(htail_id, "naca 0012.dat", AIRFOIL_DIR)
xsec_surf_id_htp = vsp.GetXSecSurf(htail_id, 0)

# root xsec keeps its own baseline value; only the outboard xsec was
# actually varied by the mesh-convergence density sweep (see main wing note above)
xsec_id_1_htp = vsp.GetXSec(xsec_surf_id_htp, 0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_1_htp, "SectTess_U"), 24.0)

xsec_id_2_htp = vsp.GetXSec(xsec_surf_id_htp, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_id_2_htp, "SectTess_U"), 20.0)

vsp.SetParmVal(htail_id, "Density", "Mass_Props", 40)
xsec_1_htp = vsp.GetXSec(xsec_surf_id_htp, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1_htp, "ThickChord"), 0.08)


# ============================================================ VERTICAL TAIL (VTP, 2 sections) ============================================================
vtp_id = vsp.AddGeom("WING")
vsp.SetGeomName(vtp_id, "VTP")
vsp.SetParmVal(vtp_id, "Tess_W", "Shape", 37)  # d=50 mesh-convergence result

vsp.InsertXSec(vtp_id, 1, vsp.XS_FOUR_SERIES)
vsp.Update()

# section 1 (root)
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

# section 2 (tip)
vsp.SetDriverGroup(vtp_id, 2, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
vsp.Update()  # required before changing section 2's drivers

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

assign_airfoil_to_component(vtp_id, "naca 0008.dat", AIRFOIL_DIR)
xsec_surf_id_vtp = vsp.GetXSecSurf(vtp_id, 0)

# root (dorsal-fin curve start) -- not varied by the mesh-convergence sweep,
# keeps its own baseline value (see main wing note above)
xsec_0_vtp = vsp.GetXSec(xsec_surf_id_vtp, 0)
vsp.SetParmVal(vsp.GetXSecParm(xsec_0_vtp, "ThickChord"), 0.08)
vsp.SetParmVal(vsp.GetXSecParm(xsec_0_vtp, "SectTess_U"), 26.0)

# junction between the curved root and the main tail section -- d=50 result
xsec_1_vtp = vsp.GetXSec(xsec_surf_id_vtp, 1)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1_vtp, "ThickChord"), 0.12)
vsp.SetParmVal(vsp.GetXSecParm(xsec_1_vtp, "SectTess_U"), 12.0)

# tail tip -- d=50 result
xsec_2_vtp = vsp.GetXSec(xsec_surf_id_vtp, 2)
vsp.SetParmVal(vsp.GetXSecParm(xsec_2_vtp, "ThickChord"), 0.08)
vsp.SetParmVal(vsp.GetXSecParm(xsec_2_vtp, "SectTess_U"), 12.0)

vsp.Update()
vsp.SetParmVal(vtp_id, "Density", "Mass_Props", 40)


# ============================================================ close open tips ============================================================
# 1 = Flat, 2 = Round
vsp.SetParmVal(wing_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(wing_id, "CapUMaxOption", "EndCap", 1.0)
vsp.SetParmVal(htail_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(htail_id, "CapUMaxOption", "EndCap", 1.0)
vsp.SetParmVal(vtp_id, "CapUMinOption", "EndCap", 1.0)
vsp.SetParmVal(vtp_id, "CapUMaxOption", "EndCap", 1.0)


# ============================================================ control surfaces (for VSPAERO) ============================================================
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


# ============================================================ VSPAERO control-surface groups ============================================================
def set_surface_gain_by_index(group_index, target_val_0, target_val_1):
    num_groups = vsp.GetNumControlSurfaceGroups()
    if group_index >= num_groups:
        print(f"Error: group {group_index} doesn't exist, there are only {num_groups} groups.")
        return

    all_settings = vsp.FindContainer("VSPAEROSettings", 0)
    all_parms = vsp.FindContainerParmIDs(all_settings)

    for pid in all_parms:
        name = vsp.GetParmName(pid)
        if "Gain" in name:
            if "_0_Gain" in name:
                vsp.SetParmVal(pid, target_val_0)
            elif "_1_Gain" in name:
                vsp.SetParmVal(pid, target_val_1)


# ailerons (differential deflection: left up, right down)
g_ailerons_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Ailerons_Group", g_ailerons_index)
vsp.AddSelectedToCSGroup([1, 2], g_ailerons_index)
set_surface_gain_by_index(0, 1, -1)

# flaps (symmetric deflection)
g_flaps_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Flaps_Group", g_flaps_index)
vsp.AddSelectedToCSGroup([3, 4], g_flaps_index)
set_surface_gain_by_index(1, 1.0, 1.0)

# elevator (symmetric deflection)
g_elevators_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Elevators_Group", g_elevators_index)
vsp.AddSelectedToCSGroup([1, 2], g_elevators_index)
set_surface_gain_by_index(2, 1.0, 1.0)

# rudder
g_rudder_index = vsp.CreateVSPAEROControlSurfaceGroup()
vsp.SetVSPAEROControlGroupName("Rudder_Group", g_rudder_index)
vsp.AddSelectedToCSGroup([3], g_rudder_index)


# ============================================================ save ============================================================
out_path = os.path.join(HERE, "AC_full_geometry_clean.vsp3")
vsp.WriteVSPFile(out_path)
print(f"plane geometry created: {out_path}")
