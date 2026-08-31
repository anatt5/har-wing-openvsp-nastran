"""
ASPECT RATIO SWEEP -- structural side, 2-SPAR CONFIGURATION.

Identical to ../AR_sweep/ar_sweep_structural.py except for the internal
FEA structure: TWO spars at 25%/75% chord (no 50% spar), matching the
final structural layout adopted in Section 7.1 (see also
../two_spar_dense_ribs_check/). The original ../AR_sweep/ run used three
spars, which was an inconsistency flagged in the thesis review -- this
folder is the corrected, 2-spar re-run of the same AR sweep.

  - Wing built with the AREA + AR + TAPER driver group (area fixed at
    0.8 m^2 total / 0.4 m^2 per panel), identical geometry convention to
    the aero sweep, for consistency.
  - Internal structure: 2 spars at 25%/75% chord, ribs spanning 5%-95% of
    the span at 9% relative spacing.
  - Skin: TU-Flex E-Glass Silenka, thickness tapered 1.0mm (root) to
    0.1mm (tip) -- held FIXED across the AR sweep (only AR itself is
    varied here; see ../../taper_sweep/ for the constant-weight
    thickness sweep at fixed AR). Spar/rib elements are unified onto the
    same Silenka material (see the "material unification" comment below),
    at their own fixed 2.5mm thickness.
  - Each AR point gets its own subfolder (AR_<value>/) with its .bdf,
    Nastran outputs, and an interactive mode-shape viewer.

Run:  python ar_sweep_structural.py
(runs MSC Nastran directly for each point -- takes a few minutes)
"""
import os
import re
import sys
import numpy as np
from scipy.spatial import cKDTree
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.op2.op2 import OP2
import plotly.graph_objects as go

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
AIRFOIL_DIR = os.path.join(REPO_ROOT, "airfoils")
NASTRAN_EXE = r"C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\msc20261\win64i8\nastran.exe"

sys.path.insert(0, REPO_ROOT)
from utils.airfoil_utils import AIRFOIL_TYPE, read_airfoil_points

# --------------------------------------------------- TU-Flex skin properties
SKIN_E1  = 36.75e9
SKIN_E2  = 12.25e9
SKIN_G12 = 3.95e9
SKIN_NU12 = 0.23
SKIN_RHO = 2550.0
T_ROOT = 1.0e-3   # m
T_TIP  = 0.1e-3   # m -- real TU-Flex value, fixed across this sweep
N_MODES = 15
S_FIXED = 0.8      # total wing area (m^2), held fixed -- matches ar_sweep_aero.py
TAPER_MW = 0.66667
SWEEP_ANGLE = 28.0


# ============================== geometry + FEA structure build ==============================
def build_structural_geometry(AR_mw, out_bdf_name):
    import openvsp as vsp
    vsp.VSPRenew()

    wing_id = vsp.AddGeom("WING", "")
    vsp.SetGeomName(wing_id, "main wing")
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", 25)

    vsp.SetDriverGroup(wing_id, 1,
                        vsp.AREA_WSECT_DRIVER,
                        vsp.AR_WSECT_DRIVER,
                        vsp.TAPER_WSECT_DRIVER)
    vsp.Update()

    vsp.SetParmVal(vsp.GetParm(wing_id, "Area", "XSec_1"), S_FIXED / 2.0)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Aspect", "XSec_1"), AR_mw)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Taper", "XSec_1"), TAPER_MW)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Sweep", "XSec_1"), SWEEP_ANGLE)
    vsp.Update()

    vsp.SetParmVal(vsp.GetParm(wing_id, "X_Rel_Location", "XForm"), 0.7)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Location", "XForm"), 0.0)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Z_Rel_Location", "XForm"), 0.0)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Y_Rel_Rotation", "XForm"), 0.0)
    vsp.Update()

    airfoil_type = AIRFOIL_TYPE
    pts_up, pts_low = read_airfoil_points(AIRFOIL_DIR, "clark Y.dat")
    xsec_surf_id_mwing = vsp.GetXSecSurf(wing_id, 0)
    for i in range(vsp.GetNumXSec(xsec_surf_id_mwing)):
        vsp.ChangeXSecShape(xsec_surf_id_mwing, i, airfoil_type)
        xsec_id = vsp.GetXSec(xsec_surf_id_mwing, i)
        vsp.SetXSecAlias(xsec_id, "clark Y")
        vsp.SetXSecCurveAlias(xsec_id, "clark Y")
        if hasattr(vsp, "vec3dVec"):
            v_up, v_low = vsp.vec3dVec(), vsp.vec3dVec()
            for p in pts_up: v_up.push_back(p)
            for p in pts_low: v_low.push_back(p)
            vsp.SetAirfoilPnts(xsec_id, v_up, v_low)
        else:
            vsp.SetAirfoilPnts(xsec_id, pts_up, pts_low)

    for i in [0, 1]:
        xsec_id = vsp.GetXSec(xsec_surf_id_mwing, i)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "SectTess_U"), 81.0)
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "Dihedral"), 6.0)

    vsp.SetParmVal(wing_id, "Density", "Mass_Props", 40)
    vsp.Update()

    # ---- internal FEA structure: 2 SPARS (25%/75%), rib array ----
    vsp.AddFeaStruct(wing_id, True)
    spar0_id = vsp.AddFeaPart(wing_id, 0, vsp.FEA_SPAR)
    rib_id   = vsp.AddFeaPart(wing_id, 0, vsp.FEA_RIB_ARRAY)
    vsp.Update()
    vsp.SetParmVal(vsp.GetParm(spar0_id, "RelCenterLocation", "FeaPart"), 0.25)

    spar1_id = vsp.AddFeaPart(wing_id, 0, vsp.FEA_SPAR)
    vsp.Update()
    vsp.SetParmVal(vsp.GetParm(spar1_id, "RelCenterLocation", "FeaPart"), 0.75)

    vsp.SetParmVal(vsp.GetParm(rib_id, "RelStartLocation", "FeaRibArray"), 0.05)
    vsp.SetParmVal(vsp.GetParm(rib_id, "RelEndLocation", "FeaRibArray"), 0.95)
    vsp.SetParmVal(vsp.GetParm(rib_id, "RibRelSpacing", "FeaRibArray"), 0.09)
    vsp.Update()

    # ---- units ----
    for c_id in vsp.FindContainers():
        for p_id in vsp.FindContainerParmIDs(c_id):
            p_name = vsp.GetParmName(p_id)
            if p_name == "StructUnit":
                vsp.SetParmVal(p_id, 0.0)
            elif p_name == "StructModelUnit":
                vsp.SetParmVal(p_id, 2.0)
    vsp.Update()

    # ---- default materials (spars/ribs keep this; skin overridden downstream) ----
    for c_id in vsp.FindContainers():
        c_name = vsp.GetContainerName(c_id)
        p_ids = vsp.FindContainerParmIDs(c_id)
        if "DefaultShell" in c_name:
            for p_id in p_ids:
                p_name = vsp.GetParmName(p_id)
                if "Mat" in p_name:
                    vsp.SetParmVal(p_id, 12.0)
                elif p_name == "Thickness":
                    vsp.SetParmVal(p_id, 0.0025)
                elif p_name == "LengthUnit":
                    vsp.SetParmVal(p_id, 2.0)
        elif "DefaultBeam" in c_name:
            for p_id in p_ids:
                p_name = vsp.GetParmName(p_id)
                if "Mat" in p_name:
                    vsp.SetParmVal(p_id, 0.0)
                elif p_name == "LengthUnit":
                    vsp.SetParmVal(p_id, 2.0)
    vsp.Update()

    # ---- FEA mesh density + export ----
    struct_ind = 0
    vsp.SetFeaMeshStructIndex(struct_ind)
    mesh_density = {"HighOrderElementFlag": 0.0, "BaseLen": 0.5, "MinLen": 0.025, "MaxGap": 0.005}
    for c_id in vsp.FindContainers():
        if vsp.GetContainerName(c_id) == "main wing_Struct0":
            for p_id in vsp.FindContainerParmIDs(c_id):
                p_name = vsp.GetParmName(p_id)
                if p_name in mesh_density:
                    vsp.SetParmVal(p_id, mesh_density[p_name])
    vsp.Update()

    vsp.SetFeaMeshFileName(wing_id, struct_ind, 1, out_bdf_name)
    vsp.Update()
    struct_id_final = vsp.GetFeaStructID(wing_id, struct_ind)
    vsp.ComputeFeaMesh(struct_id_final, 1)

    vsp3_path = os.path.splitext(out_bdf_name)[0] + ".vsp3"
    vsp.WriteVSPFile(vsp3_path)
    print(f"  written: {vsp3_path}")

    S = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalArea", "WingGeom"))
    bref = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalSpan", "WingGeom"))
    mac = vsp.GetParmVal(vsp.GetParm(wing_id, "TotalChord", "WingGeom"))
    return S, bref, mac


# ============================== mesh repair (unchanged from setup_sol103_half_clamped_tuflex.py) ==============================
def interior_angles(pts):
    pts = [np.asarray(p, dtype=float) for p in pts]
    n = len(pts)
    normal = np.zeros(3)
    for i in range(n):
        normal += np.cross(pts[i], pts[(i + 1) % n])
    norm_len = np.linalg.norm(normal)
    normal = normal / norm_len if norm_len > 1e-12 else np.array([0.0, 0.0, 1.0])
    angles = []
    for i in range(n):
        a, b, c = pts[i - 1], pts[i], pts[(i + 1) % n]
        v1, v2 = a - b, c - b
        cosang = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)
        angle = np.degrees(np.arccos(cosang))
        if np.dot(np.cross(v1, v2), normal) > 0:
            angle = 360.0 - angle
        angles.append(angle)
    return angles


def apply_node_merges(model, merges, min_tri_angle=2.0):
    def resolve(nid):
        while nid in merges:
            nid = merges[nid]
        return nid

    n_converted, n_dropped = 0, 0
    for eid, elem in list(model.elements.items()):
        if elem.type != 'CQUAD4':
            continue
        orig = elem.node_ids
        nids = [resolve(n) for n in orig]
        if nids == orig:
            continue
        uniq = []
        for n in nids:
            if not uniq or uniq[-1] != n:
                uniq.append(n)
        if len(uniq) > 1 and uniq[0] == uniq[-1]:
            uniq.pop()
        uniq = list(dict.fromkeys(uniq))
        pid = elem.pid
        del model.elements[eid]
        if len(uniq) == 4:
            model.add_cquad4(eid, pid, uniq)
        elif len(uniq) == 3:
            pts = [np.array(model.nodes[n].xyz) for n in uniq]
            if min(interior_angles(pts)) > min_tri_angle:
                model.add_ctria3(eid, pid, uniq)
                n_converted += 1
            else:
                n_dropped += 1
        else:
            n_dropped += 1

    used = set()
    for elem in model.elements.values():
        used.update(elem.node_ids)
    for nid in list(model.nodes.keys()):
        if nid not in used:
            del model.nodes[nid]
    return n_converted, n_dropped


def repair_degenerate_quads(model, angle_thresh=179.0, min_tri_angle=2.0):
    merges = {}

    def resolve(nid):
        while nid in merges:
            nid = merges[nid]
        return nid

    for _ in range(10):
        n_flagged = 0
        for elem in list(model.elements.values()):
            if elem.type != 'CQUAD4':
                continue
            nids = [resolve(n) for n in elem.node_ids]
            if len(set(nids)) < 4:
                continue
            pts = [np.array(model.nodes[n].xyz) for n in nids]
            if max(interior_angles(pts)) <= angle_thresh:
                continue
            n_flagged += 1
            closest = min(
                ((np.linalg.norm(pts[i] - pts[j]), nids[i], nids[j])
                 for i in range(4) for j in range(i + 1, 4)),
                key=lambda t: t[0],
            )
            lo, hi = sorted(closest[1:])
            if hi != lo:
                merges[hi] = lo
        if n_flagged == 0:
            break

    n_converted, n_dropped = apply_node_merges(model, merges, min_tri_angle)
    print(f"  mesh repair: merged {len(merges)} near-duplicate node pairs, "
          f"kept {n_converted} as CTRIA3, dropped {n_dropped} zero-area slivers")
    return len(merges)


def merge_weak_connectivity_nodes(model, max_conn=2, max_iter=15):
    def resolve(nid, merges):
        while nid in merges:
            nid = merges[nid]
        return nid

    merges = {}
    total_merged = 0
    for _ in range(max_iter):
        count = {}
        for elem in model.elements.values():
            for n in set(resolve(x, merges) for x in elem.node_ids):
                count[n] = count.get(n, 0) + 1
        weak = [n for n, c in count.items() if c <= max_conn]
        if not weak:
            break

        pass_merges = {}
        for n in weak:
            rn = resolve(n, merges)
            if rn in pass_merges:
                continue
            p = np.array(model.nodes[rn].xyz)
            best = None
            for elem in model.elements.values():
                ids = [resolve(x, merges) for x in elem.node_ids]
                if rn not in ids:
                    continue
                for other in ids:
                    if other == rn:
                        continue
                    d = np.linalg.norm(np.array(model.nodes[other].xyz) - p)
                    score = (count.get(other, 0), -d)
                    if best is None or score > best[0]:
                        best = (score, other)
            if best:
                pass_merges[rn] = best[1]

        if not pass_merges:
            break
        merges.update(pass_merges)
        apply_node_merges(model, merges)
        total_merged += len(pass_merges)

    print(f"  weak-connectivity repair: merged {total_merged} nodes with <= {max_conn} elements")
    return total_merged


def rigidize_remaining_weak_nodes(model, max_conn=2, root_y_tol=1e-4):
    count = {}
    for elem in model.elements.values():
        for n in elem.node_ids:
            count[n] = count.get(n, 0) + 1

    weak = [n for n, c in count.items()
            if c <= max_conn and abs(model.nodes[n].xyz[1]) > root_y_tol]

    next_eid = max(model.elements.keys(), default=0) + 1
    n_rigid = 0
    for n in weak:
        p = np.array(model.nodes[n].xyz)
        best = None
        for elem in model.elements.values():
            if n not in elem.node_ids:
                continue
            for other in elem.node_ids:
                if other == n:
                    continue
                d = np.linalg.norm(np.array(model.nodes[other].xyz) - p)
                score = (count.get(other, 0), -d)
                if best is None or score > best[0]:
                    best = (score, other)
        if best:
            model.add_rbe2(next_eid, best[1], '123456', [n])
            next_eid += 1
            n_rigid += 1

    print(f"  weak-node rigidization: added {n_rigid} RBE2 constraints")
    return n_rigid


def fix_preexisting_duplicate_id_quads(model, min_tri_angle=2.0):
    """OpenVSP's raw mesh can itself contain a CQUAD4 whose 4 corner node
    IDs are not all distinct (e.g. [977,978,979,978]) -- a degenerate
    element from the mesher, present before any of our own node merges.
    repair_degenerate_quads()/apply_node_merges() only rebuild elements
    *touched by a merge*, so a pre-existing duplicate like this slips
    through both untouched. Fix it here: convert to CTRIA3 if the 3
    distinct corners form a reasonable triangle, else drop it."""
    n_converted, n_dropped = 0, 0
    for eid, elem in list(model.elements.items()):
        if elem.type != 'CQUAD4':
            continue
        nids = elem.node_ids
        if len(set(nids)) == 4:
            continue
        uniq = list(dict.fromkeys(nids))
        pid = elem.pid
        del model.elements[eid]
        if len(uniq) == 3:
            pts = [np.array(model.nodes[n].xyz) for n in uniq]
            if min(interior_angles(pts)) > min_tri_angle:
                model.add_ctria3(eid, pid, uniq)
                n_converted += 1
            else:
                n_dropped += 1
        else:
            n_dropped += 1
    if n_converted or n_dropped:
        print(f"  pre-existing duplicate-ID quads: kept {n_converted} as CTRIA3, dropped {n_dropped}")
    return n_converted, n_dropped


def equivalence_duplicate_nodes(model, tol=2.0e-3):
    nids = list(model.nodes.keys())
    xyz = np.array([model.nodes[n].xyz for n in nids])
    pairs = cKDTree(xyz).query_pairs(tol)

    parent = {n: n for n in nids}

    def find(n):
        while parent[n] != n:
            n = parent[n]
        return n

    for i, j in pairs:
        a, b = find(nids[i]), find(nids[j])
        if a != b:
            parent[max(a, b)] = min(a, b)

    merges = {n: find(n) for n in nids if find(n) != n}
    apply_node_merges(model, merges)
    print(f"  node equivalencing: merged {len(merges)} near-duplicate nodes (tol={tol * 1e3:.1f}mm)")
    return len(merges)


def parse_skin_element_ids(comments_bdf):
    skin_ids = set()
    is_skin_block = False
    with open(comments_bdf) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('$'):
                is_skin_block = 'ShellElements' in stripped and 'Skin' in stripped
                continue
            if not stripped:
                is_skin_block = False
                continue
            if is_skin_block:
                nums = re.findall(r'\d+', stripped)
                skin_ids.update(int(n) for n in nums)
    return skin_ids


# ============================== full pipeline for one AR point ==============================
def process_one_ar(AR_mw, AR_std, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    raw_bdf = os.path.join(out_dir, "AC_wing.bdf")
    raw_bdf_bdf = raw_bdf + ".bdf"

    S, bref, mac = build_structural_geometry(AR_mw, raw_bdf)
    print(f"  S={S:.4f} m^2  bref={bref:.4f} m  mac={mac:.4f} m")

    model = BDF(debug=False)
    # validate=False: at high AR, OpenVSP's raw mesh can contain a CQUAD4
    # with two duplicate corner node IDs (e.g. [977,978,979,978]) --
    # pyNastran's own read-time validation rejects that outright before we
    # get a chance to run our own repair (which only runs after a
    # successful read). Skip that eager check here; the file is still
    # fully validated at the very end when the final, repaired bdf is
    # re-read with xref=True (default validate=True) for post-processing.
    model.read_bdf(raw_bdf_bdf, xref=True, validate=False)

    fix_preexisting_duplicate_id_quads(model)
    equivalence_duplicate_nodes(model, tol=2.0e-3)
    for _ in range(6):
        n_angle = repair_degenerate_quads(model)
        n_weak = merge_weak_connectivity_nodes(model)
        if n_angle == 0 and n_weak == 0:
            break

    skin_eids_all = parse_skin_element_ids(raw_bdf)

    y_tol = 1e-4
    drop_eids = [eid for eid, elem in model.elements.items()
                 if any(model.nodes[nid].xyz[1] > y_tol for nid in elem.node_ids)]
    for eid in drop_eids:
        del model.elements[eid]
    used = set()
    for elem in model.elements.values():
        used.update(elem.node_ids)
    for nid in list(model.nodes.keys()):
        if nid not in used:
            del model.nodes[nid]

    nids = list(model.nodes.keys())
    xyz = np.array([model.nodes[nid].xyz for nid in nids])
    print(f"  half-wing: {len(nids)} nodes | {len(model.elements)} elements")

    skin_eids = skin_eids_all & set(model.elements.keys())

    NEW_MID = max(model.materials.keys(), default=0) + 1
    model.add_mat8(NEW_MID, SKIN_E1, SKIN_E2, SKIN_NU12, g12=SKIN_G12,
                    g1z=SKIN_G12, g2z=SKIN_G12, rho=SKIN_RHO)

    y_max = np.abs(xyz[:, 1]).max()
    next_pid = max(model.properties.keys(), default=0) + 1
    for eid in skin_eids:
        elem = model.elements[eid]
        y_c = abs(np.mean([model.nodes[nid].xyz[1] for nid in elem.node_ids]))
        t = T_ROOT + (T_TIP - T_ROOT) * (y_c / y_max)
        model.add_pshell(next_pid, mid1=NEW_MID, t=t, mid2=NEW_MID, mid3=NEW_MID)
        elem.pid = next_pid
        next_pid += 1
    print(f"  skin: {len(skin_eids)} elements, taper {T_ROOT*1e3:.2f}mm->{T_TIP*1e3:.2f}mm")

    # ---- material unification: spar/rib elements also get the Silenka
    # MAT8 (same as the skin), instead of OpenVSP's uncited default
    # "Glass Epoxy" library material. Thickness (2.5mm) is untouched.
    non_skin_pids = {elem.pid for eid, elem in model.elements.items()
                      if eid not in skin_eids and elem.type in ('CQUAD4', 'CTRIA3')}
    for pid in non_skin_pids:
        prop = model.properties[pid]
        prop.mid1 = NEW_MID
        prop.mid2 = NEW_MID
        prop.mid3 = NEW_MID

    rigidize_remaining_weak_nodes(model)
    model.cross_reference()
    root_nodes = sorted(int(nid) for nid in model.nodes
                         if abs(model.nodes[nid].xyz[1]) <= y_tol)
    print(f"  root plane -> {len(root_nodes)} clamped nodes")

    model.sol = 103
    model.case_control_deck = CaseControlDeck([
        'ECHO = NONE', 'DISPLACEMENT(PLOT) = ALL', 'SPC = 1', 'METHOD = 1',
    ])
    model.spcs = {}; model.spcadds = {}; model.loads = {}; model.load_combinations = {}
    model.add_spc1(conid=1, components='123456', nodes=root_nodes)
    model.add_eigrl(sid=1, v1=None, v2=None, nd=N_MODES)
    model.add_param('POST', [-1])

    out_bdf = os.path.join(out_dir, "AC_wing_half_clamped_sol103.bdf")
    model.write_bdf(out_bdf)

    # ---- run Nastran directly ----
    import subprocess
    cmd = [NASTRAN_EXE, f"JID={out_bdf}", "mem=max"]
    subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, timeout=120)

    log_path = os.path.join(out_dir, "ac_wing_half_clamped_sol103.log")
    op2_path = os.path.join(out_dir, "ac_wing_half_clamped_sol103.op2")
    if not os.path.exists(log_path) or "NSEXIT: EXIT(0)" not in open(log_path).read():
        print(f"  [!] Nastran did not exit cleanly for AR={AR_std:.2f} -- check {log_path}")
        return None

    # ---- extract mode 1 + contamination ratio + skin mass ----
    bdf_final = BDF(debug=False)
    bdf_final.read_bdf(out_bdf, xref=True)
    op2 = OP2(debug=False)
    op2.read_op2(op2_path)
    eig = op2.eigenvectors[1]
    freqs_hz = eig.mode_cycles
    disp = eig.data[0, :, :3]
    mags = np.linalg.norm(disp, axis=1)
    med = np.median(mags[mags > 1e-12])
    ratio = mags.max() / med if med > 0 else float('inf')

    skin_mass = other_mass = 0.0
    for elem in bdf_final.elements.values():
        if elem.type not in ('CQUAD4', 'CTRIA3'):
            continue
        if elem.eid in skin_eids:
            skin_mass += elem.Mass()
        else:
            other_mass += elem.Mass()

    # interactive mode-shape viewer for this AR point
    write_interactive_html(out_dir, bdf_final, eig, skin_eids)

    return dict(AR_mw=AR_mw, AR_std=AR_std, S=S, bref=bref, mac=mac,
                f1_hz=freqs_hz[0], f1_ratio=ratio,
                skin_mass_g=skin_mass * 1000, other_mass_g=other_mass * 1000,
                total_mass_g=(skin_mass + other_mass) * 1000)


def write_interactive_html(folder, bdf_model, eig, skin_eids, out_name="mode_shapes_interactive.html"):
    node_ids = eig.node_gridtype[:, 0]
    node_index = {nid: i for i, nid in enumerate(node_ids)}
    xyz0 = np.array([bdf_model.nodes[nid].xyz for nid in node_ids])
    faces = [elem.node_ids for elem in bdf_model.elements.values()
             if elem.type in ('CQUAD4', 'CTRIA3') and elem.eid in skin_eids]
    faces_idx = [[node_index[n] for n in f] for f in faces]

    tri_i, tri_j, tri_k = [], [], []
    edges_i, edges_j = [], []
    for f in faces_idx:
        if len(f) == 4:
            tri_i += [f[0], f[0]]; tri_j += [f[1], f[2]]; tri_k += [f[2], f[3]]
        else:
            tri_i.append(f[0]); tri_j.append(f[1]); tri_k.append(f[2])
        for a, b in zip(f, f[1:] + f[:1]):
            edges_i.append(a); edges_j.append(b)

    n_modes = eig.data.shape[0]
    freqs_hz = eig.mode_cycles
    span = xyz0[:, 1].max() - xyz0[:, 1].min()
    DISP_TARGET_FRAC = 0.25

    wf_x, wf_y, wf_z = [], [], []
    for a, b in zip(edges_i, edges_j):
        wf_x += [xyz0[a, 0], xyz0[b, 0], None]
        wf_y += [xyz0[a, 1], xyz0[b, 1], None]
        wf_z += [xyz0[a, 2], xyz0[b, 2], None]

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=wf_x, y=wf_y, z=wf_z, mode='lines',
                                line=dict(color='#9aa5b1', width=2.5), opacity=0.6,
                                name='undeformed', hoverinfo='skip', showlegend=False))
    for m in range(n_modes):
        d = eig.data[m, :, :3]
        mags = np.linalg.norm(d, axis=1)
        max_disp = mags.max()
        scale = DISP_TARGET_FRAC * span / max_disp if max_disp > 0 else 1.0
        xyz_def = xyz0 + scale * d
        intensity = mags / max_disp if max_disp > 0 else mags
        fig.add_trace(go.Mesh3d(
            x=xyz_def[:, 0], y=xyz_def[:, 1], z=xyz_def[:, 2], i=tri_i, j=tri_j, k=tri_k,
            intensity=intensity, intensitymode='vertex', colorscale='Turbo', cmin=0, cmax=1,
            colorbar=dict(title='Rel. disp.', thickness=18, len=0.75, x=0.83, xanchor='left',
                           y=0.5, yanchor='middle', xpad=0) if m == 0 else None,
            showscale=(m == 0), opacity=0.98, flatshading=False,
            name=f'Mode {m + 1}', visible=(m == 0),
            lighting=dict(ambient=0.55, diffuse=0.85, specular=0.3, roughness=0.5),
            lightposition=dict(x=100, y=200, z=150)))

    buttons = [dict(label=f"Mode {m + 1}: {freqs_hz[m]:.2f} Hz", method="update",
                     args=[{"visible": [True] + [i == m for i in range(n_modes)]},
                           {"title": f"Mode {m + 1}: {freqs_hz[m]:.2f} Hz"}])
               for m in range(n_modes)]

    fig.update_layout(
        template='plotly_white',
        title=dict(text=f"Mode 1: {freqs_hz[0]:.2f} Hz", font=dict(size=20, color='#2c3e50')),
        scene=dict(aspectmode='data', domain=dict(x=[0, 0.8], y=[0, 1]),
                   xaxis=dict(title='X', backgroundcolor='#f7f9fb', gridcolor='#dfe6ec'),
                   yaxis=dict(title='Y (span)', backgroundcolor='#f7f9fb', gridcolor='#dfe6ec'),
                   zaxis=dict(title='Z', backgroundcolor='#f7f9fb', gridcolor='#dfe6ec')),
        paper_bgcolor='#ffffff',
        updatemenus=[dict(buttons=buttons, direction="down", showactive=True,
                           x=0.02, y=0.98, xanchor='left', yanchor='top',
                           bgcolor='#eef2f6', bordercolor='#c3cdd6')],
        margin=dict(l=0, r=0, t=50, b=0), height=850)
    fig.write_html(os.path.join(folder, out_name))


# ============================== MAIN ==============================
if __name__ == "__main__":
    AR_points = [(3.2, 6.4), (4.8, 9.6), (6.4, 12.8), (8.0, 16.0), (9.6, 19.2)]
    results = []
    for AR_mw, AR_std in AR_points:
        print(f"\n{'#'*70}\n#  AR_mw={AR_mw:.2f}  AR_std={AR_std:.2f}\n{'#'*70}")
        out_dir = os.path.join(HERE, f"AR_{AR_std:.2f}")
        r = process_one_ar(AR_mw, AR_std, out_dir)
        if r is not None:
            print(f"  f1={r['f1_hz']:.3f} Hz  ratio={r['f1_ratio']:.2f}  "
                  f"skin_mass={r['skin_mass_g']:.1f}g  other_mass={r['other_mass_g']:.1f}g")
            results.append(r)

    out_csv = os.path.join(HERE, "ar_sweep_structural_results.csv")
    with open(out_csv, "w") as f:
        f.write("AR_mw,AR_standard,S,bref,mac,f1_hz,f1_ratio,skin_mass_g,other_mass_g,total_mass_g\n")
        for r in results:
            f.write(f"{r['AR_mw']:.3f},{r['AR_std']:.3f},{r['S']:.4f},{r['bref']:.4f},{r['mac']:.4f},"
                    f"{r['f1_hz']:.4f},{r['f1_ratio']:.2f},{r['skin_mass_g']:.2f},"
                    f"{r['other_mass_g']:.2f},{r['total_mass_g']:.2f}\n")
    print(f"\nwritten: {out_csv}")
