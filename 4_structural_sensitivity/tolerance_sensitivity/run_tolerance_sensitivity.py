"""
NODE-EQUIVALENCING TOLERANCE SENSITIVITY -- addresses the thesis review
comment: better justify the tolerance used to merge near-duplicate
nodes (equivalence_duplicate_nodes(), Section~\\ref{sec: wing_structure}).

Reuses the SAME raw mesh already built for the BaseLen=0.5 case of
../mesh_convergence/ (identical 2-spar geometry, no need to rebuild it in
OpenVSP again) and re-runs the full repair + TU-Flex taper + half-wing
clamp + SOL 103 pipeline from scratch for each candidate tolerance,
tracking:
  - how many node pairs the equivalencing step itself merges at that
    tolerance (too small -> some real duplicate seams are missed; too
    large -> risk of merging genuinely distinct nodes)
  - the resulting first bending frequency

Run:  python run_tolerance_sensitivity.py
"""
import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.op2.op2 import OP2
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(os.path.dirname(HERE), "mesh_convergence", "BaseLen_0.50")
RAW_BDF = os.path.join(SRC_DIR, "AC_wing.bdf")
RAW_BDF_BDF = os.path.join(SRC_DIR, "AC_wing.bdf.bdf")
NASTRAN_EXE = r"C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\msc20261\win64i8\nastran.exe"

SKIN_E1, SKIN_E2, SKIN_G12, SKIN_NU12, SKIN_RHO = 36.75e9, 12.25e9, 3.95e9, 0.23, 2550.0
T_ROOT, T_TIP = 1.0e-3, 0.1e-3
N_MODES = 5

# candidate tolerances (mm), current thesis-wide choice (2mm) among them
TOLERANCES_MM = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
THESIS_TOL_MM = 2.0


# ============================== mesh repair helpers (same as run_mesh_convergence.py) ==============================
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


# ============================== full pipeline for one tolerance value ==============================
def process_one_tolerance(tol_mm, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    model = BDF(debug=False)
    model.read_bdf(RAW_BDF_BDF, xref=True, validate=False)

    fix_preexisting_duplicate_id_quads(model)
    n_equiv_merged = equivalence_duplicate_nodes(model, tol=tol_mm * 1e-3)
    for _ in range(6):
        n_angle = repair_degenerate_quads(model)
        n_weak = merge_weak_connectivity_nodes(model)
        if n_angle == 0 and n_weak == 0:
            break

    skin_eids_all = parse_skin_element_ids(RAW_BDF)

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

    n_rigid = rigidize_remaining_weak_nodes(model)
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
        print(f"  [!] Nastran did not exit cleanly for tol={tol_mm}mm -- check {log_path}")
        return None

    op2 = OP2(debug=False)
    op2.read_op2(op2_path)
    eig = op2.eigenvectors[1]
    f1_hz = float(eig.mode_cycles[0])

    return dict(tol_mm=tol_mm, n_equiv_merged=n_equiv_merged, n_rigid=n_rigid,
                n_nodes=n_nodes, n_elements=n_elements, f1_hz=f1_hz)


# ============================== MAIN ==============================
if __name__ == "__main__":
    results = []
    for tol_mm in TOLERANCES_MM:
        label = f"tol_{tol_mm:.1f}mm"
        out_dir = os.path.join(HERE, label)
        print(f"\n{'#'*70}\n#  {label}\n{'#'*70}")
        try:
            r = process_one_tolerance(tol_mm, out_dir)
        except Exception as exc:
            # at very large tolerances (e.g. 50mm) the equivalencing step
            # collapses/degenerates the mesh so badly that downstream steps
            # fail outright (same failure mode as the original,
            # uncited-material tolerance sweep, which likewise has no
            # tol=50.0mm row in its results) -- exclude this point rather
            # than aborting the whole sweep.
            print(f"  [!] tol={tol_mm}mm crashed ({exc!r}) -- excluding this point, same as the original study")
            r = None
        if r is not None:
            print(f"  merged={r['n_equiv_merged']}  rigidized={r['n_rigid']}  "
                  f"n_elements={r['n_elements']}  f1={r['f1_hz']:.4f} Hz")
            results.append(r)

    out_csv = os.path.join(HERE, "tolerance_sensitivity_results.csv")
    with open(out_csv, "w") as f:
        f.write("tol_mm,n_equiv_merged,n_rigidized,n_nodes,n_elements,f1_hz\n")
        for r in results:
            f.write(f"{r['tol_mm']:.2f},{r['n_equiv_merged']},{r['n_rigid']},"
                    f"{r['n_nodes']},{r['n_elements']},{r['f1_hz']:.4f}\n")
    print(f"\nwritten: {out_csv}")

    if results:
        tols = [r['tol_mm'] for r in results]
        merged = [r['n_equiv_merged'] for r in results]
        f1s = [r['f1_hz'] for r in results]

        plt.rcParams.update({"font.size": 12})
        fig, ax1 = plt.subplots(figsize=(7, 5))
        ax2 = ax1.twinx()
        l1, = ax1.plot(tols, merged, "o-", color="#c0392b", linewidth=2, markersize=7, label="Nodes merged by equivalencing")
        l2, = ax2.plot(tols, f1s, "s--", color="#1f6fb2", linewidth=2, markersize=6, label="First bending frequency")
        ax1.set_xscale("log")
        ax1.set_xlabel("Node-equivalencing tolerance (mm)")
        ax1.set_ylabel("Nodes merged by this step", color="#c0392b")
        ax2.set_ylabel("First bending frequency of the clamped\nhalf-wing structural model (Hz)", color="#1f6fb2")
        ax1.tick_params(axis='y', labelcolor="#c0392b")
        ax2.tick_params(axis='y', labelcolor="#1f6fb2")
        ax1.axvline(THESIS_TOL_MM, color="#7f8c8d", linestyle=":", linewidth=1.5)
        ax1.annotate("value used\nin this thesis", (THESIS_TOL_MM, max(merged)*0.9), fontsize=8, color="#555555", ha="center")
        ax1.set_title("Sensitivity to the node-equivalencing tolerance")
        ax1.grid(alpha=0.3)
        ax1.legend(handles=[l1, l2], loc="center right", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "fig_tolerance_sensitivity.png"), dpi=160)
        plt.close(fig)
        print("written: fig_tolerance_sensitivity.png")
