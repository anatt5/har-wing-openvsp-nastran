"""
Shared OpenVSP airfoil-assignment helpers.

Every geometry-building script in this repo (the full-aircraft geometry,
the VSPAERO studies, and the structural mesh studies) needs to load an
airfoil's coordinates from a .dat file and assign them to a WING
geometry's cross-sections. This used to be copy-pasted into each script
separately; it now lives here once so all of them stay in sync.

Each script still owns its own AIRFOIL_DIR (the path to the shared
repo-root "airfoils/" folder, computed relative to that script's own
location) and passes it in explicitly -- this module has no opinion on
where the .dat files live.
"""
import os

import openvsp as vsp


def find_airfoil_type_index():
    """OpenVSP's "airfoil from file" xsec-shape constant is named
    differently across versions -- try both known names."""
    for name in ["AF_FILE", "XS_FILE_AIRFOIL"]:
        if hasattr(vsp, name):
            return getattr(vsp, name)
    return 5


AIRFOIL_TYPE = find_airfoil_type_index()


def read_airfoil_points(airfoil_dir, file_name):
    """Read an airfoil .dat file and split it into upper/lower surface points."""
    path = os.path.join(airfoil_dir, file_name)
    up, low = [], []
    reached_le = False
    try:
        with open(path, 'r') as f:
            lines = f.readlines()[1:]  # skip the airfoil-name header line
            for line in lines:
                val = line.split()
                if len(val) == 2:
                    p = vsp.vec3d(float(val[0]), float(val[1]), 0.0)
                    if not reached_le:
                        up.append(p)
                        if float(val[0]) <= 0.0:
                            reached_le = True
                    else:
                        low.append(p)
        up.reverse()
        return up, low
    except FileNotFoundError:
        print(f"Error: could not find {file_name} at {path}")
        return None, None


def assign_airfoil_to_component(geom_id, dat_name, airfoil_dir):
    """Assign an airfoil .dat file's coordinates to every cross-section of geom_id."""
    pts_up, pts_low = read_airfoil_points(airfoil_dir, dat_name)
    if pts_up is None:
        return

    alias = os.path.splitext(dat_name)[0]
    xsec_surf_id = vsp.GetXSecSurf(geom_id, 0)
    n_xsecs = vsp.GetNumXSec(xsec_surf_id)

    for i in range(n_xsecs):
        vsp.ChangeXSecShape(xsec_surf_id, i, AIRFOIL_TYPE)
        xsec_id = vsp.GetXSec(xsec_surf_id, i)
        vsp.SetXSecAlias(xsec_id, alias)
        vsp.SetXSecCurveAlias(xsec_id, alias)
        if hasattr(vsp, "vec3dVec"):
            v_up, v_low = vsp.vec3dVec(), vsp.vec3dVec()
            for p in pts_up: v_up.push_back(p)
            for p in pts_low: v_low.push_back(p)
            vsp.SetAirfoilPnts(xsec_id, v_up, v_low)
        else:
            vsp.SetAirfoilPnts(xsec_id, pts_up, pts_low)
