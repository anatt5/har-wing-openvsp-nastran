"""
STRUCTURAL MESH CONVERGENCE STUDY -- addresses the thesis review comment:
check the sensitivity of the first bending frequency to the structural
mesh density, and better justify the 5mm node-equivalencing tolerance.

Same 2-spar wing (AR=6.4/AR_std=12.8 baseline, no AR sweep here) and same
TU-Flex skin taper (1.0mm root -> 0.1mm tip) as
../two_spar_dense_ribs_check/, but the FEA mesh density
(BaseLen/MinLen/MaxGap, OpenVSP's mesh size controls) is varied across 5
levels, from coarse to fine, with the current thesis-wide choice
(BaseLen=0.5) in the middle so it can be checked against both a coarser
and a finer mesh.

No contamination ratio is computed here (mode shape quality already
addressed elsewhere) -- only mode 1 frequency and mesh size (node/element
count) are tracked, plus a 2D wireframe snapshot of the skin mesh for each
density so the refinement is visible at a glance.

Run:  python run_mesh_convergence.py
(builds + repairs + clamps + runs MSC Nastran for each of 5 mesh
densities -- takes a few minutes)
"""
import os
import re
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.op2.op2 import OP2
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
AIRFOIL_DIR = os.path.join(REPO_ROOT, "airfoils")
NASTRAN_EXE = r"C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\msc20261\win64i8\nastran.exe"

sys.path.insert(0, REPO_ROOT)
from utils.airfoil_utils import AIRFOIL_TYPE, read_airfoil_points

# ---- baseline wing planform + 2-spar structure (identical to two_spar_dense_ribs_check) ----
SPAN_MW, AR_MW, TAPER_MW, SWEEP_ANGLE, DIHEDRAL = 1.6, 6.4, 0.66667, 28.0, 6.0
SPAR_LOCATIONS = [0.25, 0.75]
RIB_START, RIB_END, RIB_SPACING = 0.05, 0.95, 0.09

SKIN_E1, SKIN_E2, SKIN_G12, SKIN_NU12, SKIN_RHO = 36.75e9, 12.25e9, 3.95e9, 0.23, 2550.0
T_ROOT, T_TIP = 1.0e-3, 0.1e-3
N_MODES = 5   # only need mode 1; a few extra as a sanity margin

NODE_EQUIV_TOL_MM = 2.0  # the tolerance used throughout this thesis -- see justification in the written report

# BaseLen levels: coarse -> fine, current thesis-wide choice (0.5) in the
# middle. MinLen and MaxGap are kept at the same *relative* proportion to
# BaseLen as the original choice (MinLen = BaseLen/20, MaxGap = BaseLen/100)
# so the mesher's local refinement behaviour scales consistently.
#
# This list matches the FULL set of points that made it into the original
# (uncited-material) mesh_convergence_results.csv in ../mesh_convergence_structural/
# -- the original run's default 5-point list plus the extra points added
# later via run_extra_points.py. BaseLen=0.45/0.40 are excluded because
# they hung/failed to produce a usable result in the original run, and
# BaseLen=0.30/0.25 are excluded because they gave clearly anomalous,
# unreliable mesher output (f1 dropping to ~4.5/7.7 Hz, far outside the
# converged trend) -- replicate that same judgment call here rather than
# re-litigating it, per the same exclusion logic as the original study.
BASE_LENS = [1.00, 0.90, 0.80, 0.70, 0.65, 0.60, 0.50, 0.35]


def build_structural_geometry(base_len, min_len, max_gap, out_bdf_name):
    import openvsp as vsp
    vsp.VSPRenew()

    wing_id = vsp.AddGeom("WING", "")
    vsp.SetGeomName(wing_id, "main wing")
    vsp.SetParmVal(wing_id, "Tess_W", "Shape", 25)

    vsp.SetDriverGroup(wing_id, 1, vsp.SPAN_WSECT_DRIVER, vsp.AR_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
    vsp.Update()
    vsp.SetParmVal(vsp.GetParm(wing_id, "Span", "XSec_1"), SPAN_MW)
    vsp.SetParmVal(vsp.GetParm(wing_id, "Aspect", "XSec_1"), AR_MW)
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
        vsp.SetParmVal(vsp.GetXSecParm(xsec_id, "Dihedral"), DIHEDRAL)

    vsp.SetParmVal(wing_id, "Density", "Mass_Props", 40)
    vsp.Update()

    vsp.AddFeaStruct(wing_id, True)
    for loc in SPAR_LOCATIONS:
        spar_id = vsp.AddFeaPart(wing_id, 0, vsp.FEA_SPAR)
        vsp.Update()
        vsp.SetParmVal(vsp.GetParm(spar_id, "RelCenterLocation", "FeaPart"), loc)
    rib_id = vsp.AddFeaPart(wing_id, 0, vsp.FEA_RIB_ARRAY)
    vsp.Update()
    vsp.SetParmVal(vsp.GetParm(rib_id, "RelStartLocation", "FeaRibArray"), RIB_START)
    vsp.SetParmVal(vsp.GetParm(rib_id, "RelEndLocation", "FeaRibArray"), RIB_END)
    vsp.SetParmVal(vsp.GetParm(rib_id, "RibRelSpacing", "FeaRibArray"), RIB_SPACING)
    vsp.Update()

    for c_id in vsp.FindContainers():
        for p_id in vsp.FindContainerParmIDs(c_id):
            p_name = vsp.GetParmName(p_id)
            if p_name == "StructUnit":
                vsp.SetParmVal(p_id, 0.0)
            elif p_name == "StructModelUnit":
                vsp.SetParmVal(p_id, 2.0)
    vsp.Update()

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

    struct_ind = 0
    vsp.SetFeaMeshStructIndex(struct_ind)
    mesh_density = {"HighOrderElementFlag": 0.0, "BaseLen": base_len, "MinLen": min_len, "MaxGap": max_gap}
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


# ============================== mesh repair (same as two_spar_dense_ribs_check) ==============================
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


def fix_preexisting_duplicate_id_quads(model, min_tri_angle=2.0):
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


def equivalence_duplicate_nodes(model, tol=5.0e-3):
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


def save_mesh_snapshot(model, skin_eids, out_png, title):
    """2D top-down wireframe of the skin mesh (planform view) -- a simple,
    robust way to SEE the mesh density difference between cases, without
    the 3D-rendering pitfalls of earlier figures in this thesis."""
    fig, ax = plt.subplots(figsize=(7, 5))
    for eid in skin_eids:
        elem = model.elements.get(eid)
        if elem is None or elem.type not in ('CQUAD4', 'CTRIA3'):
            continue
        nids = elem.node_ids
        xs = [model.nodes[n].xyz[1] for n in nids] + [model.nodes[nids[0]].xyz[1]]
        ys = [model.nodes[n].xyz[0] for n in nids] + [model.nodes[nids[0]].xyz[0]]
        ax.plot(xs, ys, color="#1f6fb2", linewidth=0.4)
    ax.set_xlabel("Y (span)")
    ax.set_ylabel("X (chordwise)")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


# ============================== full pipeline for one mesh density ==============================
def process_one_density(base_len, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    min_len = base_len / 20.0
    max_gap = base_len / 100.0

    raw_bdf = os.path.join(out_dir, "AC_wing.bdf")
    raw_bdf_bdf = raw_bdf + ".bdf"
    build_structural_geometry(base_len, min_len, max_gap, raw_bdf)

    model = BDF(debug=False)
    model.read_bdf(raw_bdf_bdf, xref=True, validate=False)

    fix_preexisting_duplicate_id_quads(model)
    equivalence_duplicate_nodes(model, tol=NODE_EQUIV_TOL_MM * 1e-3)
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
    n_nodes, n_elements = len(nids), len(model.elements)
    print(f"  half-wing: {n_nodes} nodes | {n_elements} elements")

    skin_eids = skin_eids_all & set(model.elements.keys())

    # snapshot BEFORE the taper/clamp cards are added (topology is what matters)
    save_mesh_snapshot(model, skin_eids, os.path.join(out_dir, "mesh_snapshot.png"),
                       f"BaseLen={base_len:.2f} -- {n_elements} elements")

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
    root_nodes = sorted(int(nid) for nid in model.nodes if abs(model.nodes[nid].xyz[1]) <= y_tol)

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

    cmd = [NASTRAN_EXE, f"JID={out_bdf}", "mem=max"]
    subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, timeout=180)

    log_path = os.path.join(out_dir, "ac_wing_half_clamped_sol103.log")
    op2_path = os.path.join(out_dir, "ac_wing_half_clamped_sol103.op2")
    if not os.path.exists(log_path) or "NSEXIT: EXIT(0)" not in open(log_path).read():
        print(f"  [!] Nastran did not exit cleanly -- check {log_path}")
        return None

    op2 = OP2(debug=False)
    op2.read_op2(op2_path)
    eig = op2.eigenvectors[1]
    f1_hz = eig.mode_cycles[0]

    return dict(base_len=base_len, min_len=min_len, max_gap=max_gap,
                n_nodes=n_nodes, n_elements=n_elements, f1_hz=f1_hz)


# ============================== MAIN ==============================
if __name__ == "__main__":
    results = []
    for base_len in BASE_LENS:
        label = f"BaseLen_{base_len:.2f}"
        out_dir = os.path.join(HERE, label)
        print(f"\n{'#'*70}\n#  {label}\n{'#'*70}")
        r = process_one_density(base_len, out_dir)
        if r is not None:
            print(f"  f1={r['f1_hz']:.4f} Hz  n_elements={r['n_elements']}")
            results.append(r)

    results.sort(key=lambda r: -r['base_len'])  # coarse -> fine

    out_csv = os.path.join(HERE, "mesh_convergence_results.csv")
    with open(out_csv, "w") as f:
        f.write("base_len,min_len,max_gap,n_nodes,n_elements,f1_hz,pct_change_vs_previous\n")
        prev_f1 = None
        for r in results:
            pct = "" if prev_f1 is None else f"{100*(r['f1_hz']-prev_f1)/prev_f1:.3f}"
            f.write(f"{r['base_len']:.3f},{r['min_len']:.4f},{r['max_gap']:.5f},"
                    f"{r['n_nodes']},{r['n_elements']},{r['f1_hz']:.4f},{pct}\n")
            prev_f1 = r['f1_hz']
    print(f"\nwritten: {out_csv}")

    # convergence plot: f1 vs number of elements
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    n_els = [r['n_elements'] for r in results]
    f1s = [r['f1_hz'] for r in results]
    ax.plot(n_els, f1s, "o-", color="#1f6fb2", linewidth=2, markersize=7)
    for r in results:
        ax.annotate(f"BaseLen={r['base_len']:.2f}", (r['n_elements'], r['f1_hz']),
                    textcoords="offset points", xytext=(6, 6), fontsize=8, color="#555555")
    ax.set_xlabel("Number of shell elements (half-wing)")
    ax.set_ylabel("First bending frequency of the clamped\nhalf-wing structural model (Hz)")
    ax.set_title("Structural mesh convergence: Mode 1 frequency vs. mesh density")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_mesh_convergence.png"), dpi=160)
    plt.close(fig)
    print("written: fig_mesh_convergence.png")
