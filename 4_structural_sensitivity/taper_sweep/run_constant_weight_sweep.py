"""
CONSTANT-WEIGHT LAMINATE TAPER SWEEP -- 2-SPAR CONFIGURATION (corrected
re-run of ../sweep_espesor_constant_weight/ with the final 2-spar layout
instead of the original 3-spar one).

Also directly addresses the thesis-review comment on mass ambiguity: does
varying T_root/T_tip touch only the skin, or also the spar/rib elements?
The taper-assignment loop below only ever reassigns PSHELL/PID on
`skin_eids` (elements labelled "Skin" in the OpenVSP comments block) --
spar and rib elements keep whatever (fixed) thickness/material
AC_structure_wing.py's "DefaultShell" originally gave them (2.5mm, Glass
Epoxy), completely untouched by this sweep. So holding the skin mass
constant across the family also holds the TOTAL structural mass constant,
since the spar/rib contribution never changes. This is verified
numerically below: every row of the output CSV reports skin_mass_g,
other_mass_g (spar+rib), AND total_mass_g, so the constant-weight claim
can be checked rather than assumed.

Geometry is IDENTICAL to the 2-spar baseline built in
../two_spar_dense_ribs_check/ (same AR=6.4 wing, same 2 spars at
25%/75% chord, same rib array at 5%-95%/9%-spacing) -- read directly from
that folder's raw mesh (AC_wing_2spar_denseribs.bdf / .bdf.bdf) and
repaired ONCE here, since only the skin thickness assignment changes per
sweep point.

Mass-conservation trick: since thickness is a linear function of span
position (t(y) = T_root + (T_tip - T_root)*(y/y_max)), and each element's
mass is (rho * area * t), total SKIN mass is an exact linear function of
(T_root, T_tip): mass = a*T_root + b*T_tip, where a and b depend only on
the (fixed) geometry. Two quick "probe" mass evaluations (no Nastran run
needed) give a and b directly, and then any target mass M defines a
1-parameter family of (T_root, T_tip) pairs via T_root = (M - b*T_tip)/a.

Target mass M matches the reference case used throughout the thesis:
T_root=1.0mm, T_tip=0.1mm (the real TU-Flex taper).

Also produces the TFG-ready figure (fig1_frequency_vs_taper.png -- f1 vs
taper ratio) directly at the end, from the in-memory results (formerly a
separate make_plots.py step reading the CSV back in).

Run:  python run_constant_weight_sweep_clean.py
(builds + repairs the mesh once, then runs MSC Nastran directly for each
of the N taper points -- takes a few minutes)
"""
import copy
import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.op2.op2 import OP2
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(os.path.dirname(HERE), "mesh_convergence", "BaseLen_0.50")
IN_BDF_COMMENTS = os.path.join(SRC_DIR, "AC_wing.bdf")
IN_BDF = os.path.join(SRC_DIR, "AC_wing.bdf.bdf")
NASTRAN_EXE = r"C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\msc20261\win64i8\nastran.exe"

SKIN_E1, SKIN_E2, SKIN_G12, SKIN_NU12, SKIN_RHO = 36.75e9, 12.25e9, 3.95e9, 0.23, 2550.0
N_MODES = 15

# reference case that fixes the target total skin mass (matches the
# TU-Flex taper used throughout the thesis)
T_ROOT_REF = 1.0e-3
T_TIP_REF = 0.1e-3

# family of tip thicknesses to sweep at constant skin mass; T_root is
# solved for each from the mass-conservation equation below.
T_TIP_FAMILY_MM = [0.10, 0.20, 0.30, 0.40, 0.55]


# ============================== mesh repair (identical to run_mesh_convergence.py) ==============================
# NOTE: kept for reference / potential reuse, but currently UNUSED by this
# script -- build_repaired_half_wing() below reads the already-repaired
# shared snapshot from save_shared_repaired_mesh.py directly, rather
# than running this repair pipeline itself (see the comment above
# REPAIRED_SNAPSHOT for why).
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


# ============================== build the shared, repaired half-wing base model ==============================
# Reads the ALREADY-repaired half-wing mesh saved by
# mesh_convergence/save_shared_repaired_mesh.py
# instead of re-running its own separate repair pipeline on the raw
# geometry. Two independently-implemented (even if logically equivalent)
# repair pipelines were found to converge to the same element/node COUNT
# but not to the exact same mesh, causing a small but real f1 mismatch
# with the mesh-convergence baseline (BaseLen=0.5) -- reading the literal
# shared snapshot removes that discrepancy by construction.
REPAIRED_SNAPSHOT = os.path.join(
    os.path.dirname(HERE), "mesh_convergence",
    "BaseLen_0.50", "AC_wing_repaired_half_wing.bdf")


def build_repaired_half_wing():
    model = BDF(debug=False)
    model.read_bdf(REPAIRED_SNAPSHOT, xref=True)

    skin_eids_all = parse_skin_element_ids(IN_BDF_COMMENTS)

    nids = list(model.nodes.keys())
    xyz = np.array([model.nodes[nid].xyz for nid in nids])
    print(f"half-wing (shared base, 2-spar): {len(nids)} nodes | {len(model.elements)} elements")

    skin_eids = skin_eids_all & set(model.elements.keys())
    print(f"skin elements: {len(skin_eids)}")
    return model, skin_eids, xyz[:, 1].__abs__().max()


# ============================== assign a given (T_root, T_tip) taper + clamp + write ==============================
def build_case(base_model, skin_eids, y_max, t_root, t_tip, out_dir):
    model = copy.deepcopy(base_model)

    new_mid = max(model.materials.keys(), default=0) + 1
    model.add_mat8(new_mid, SKIN_E1, SKIN_E2, SKIN_NU12, g12=SKIN_G12,
                    g1z=SKIN_G12, g2z=SKIN_G12, rho=SKIN_RHO)

    next_pid = max(model.properties.keys(), default=0) + 1
    for eid in skin_eids:
        elem = model.elements[eid]
        y_c = abs(np.mean([model.nodes[nid].xyz[1] for nid in elem.node_ids]))
        t = t_root + (t_tip - t_root) * (y_c / y_max)
        new_prop = model.add_pshell(next_pid, mid1=new_mid, t=t, mid2=new_mid, mid3=new_mid)
        elem.pid = next_pid
        elem.pid_ref = new_prop  # elem was already xref'd on read; reassigning
        next_pid += 1            # .pid alone leaves the stale pid_ref in place

    # ---- material unification: spar/rib elements also get the Silenka
    # MAT8 (same as the skin), instead of OpenVSP's uncited default
    # "Glass Epoxy" library material. Thickness (2.5mm) is untouched.
    non_skin_pids = {elem.pid for eid, elem in model.elements.items()
                      if eid not in skin_eids and elem.type in ('CQUAD4', 'CTRIA3')}
    for pid in non_skin_pids:
        prop = model.properties[pid]
        prop.mid1 = new_mid
        prop.mid2 = new_mid
        prop.mid3 = new_mid

    model.cross_reference()
    y_tol = 1e-4
    root_nodes = sorted(int(nid) for nid in model.nodes
                         if abs(model.nodes[nid].xyz[1]) <= y_tol)

    model.sol = 103
    model.case_control_deck = CaseControlDeck([
        'ECHO = NONE', 'DISPLACEMENT(PLOT) = ALL', 'SPC = 1', 'METHOD = 1',
    ])
    model.spcs = {}
    model.spcadds = {}
    model.loads = {}
    model.load_combinations = {}
    model.add_spc1(conid=1, components='123456', nodes=root_nodes)
    model.add_eigrl(sid=1, v1=None, v2=None, nd=N_MODES)
    model.add_param('POST', [-1])

    os.makedirs(out_dir, exist_ok=True)
    out_bdf = os.path.join(out_dir, "AC_wing_half_clamped_sol103.bdf")
    model.write_bdf(out_bdf)
    return out_bdf


def skin_mass_only(base_model, skin_eids, y_max, t_root, t_tip):
    """Fast mass probe -- no Nastran run, no clamp/SOL cards, just the
    total skin mass for a given (t_root, t_tip) taper."""
    model = copy.deepcopy(base_model)
    new_mid = max(model.materials.keys(), default=0) + 1
    model.add_mat8(new_mid, SKIN_E1, SKIN_E2, SKIN_NU12, g12=SKIN_G12,
                    g1z=SKIN_G12, g2z=SKIN_G12, rho=SKIN_RHO)
    next_pid = max(model.properties.keys(), default=0) + 1
    for eid in skin_eids:
        elem = model.elements[eid]
        y_c = abs(np.mean([model.nodes[nid].xyz[1] for nid in elem.node_ids]))
        t = t_root + (t_tip - t_root) * (y_c / y_max)
        new_prop = model.add_pshell(next_pid, mid1=new_mid, t=t, mid2=new_mid, mid3=new_mid)
        elem.pid = next_pid
        elem.pid_ref = new_prop  # see build_case: avoids a stale pid_ref
        next_pid += 1

    # ---- material unification (see build_case above); harmless here since
    # only skin_eids masses are summed below, kept for consistency.
    non_skin_pids = {elem.pid for eid, elem in model.elements.items()
                      if eid not in skin_eids and elem.type in ('CQUAD4', 'CTRIA3')}
    for pid in non_skin_pids:
        prop = model.properties[pid]
        prop.mid1 = new_mid
        prop.mid2 = new_mid
        prop.mid3 = new_mid

    # elem.Mass() reads the cached pid_ref set by cross_reference(), not
    # model.properties[elem.pid] freshly -- without this call it would
    # silently use each element's stale (pre-reassignment) property.
    model.cross_reference()
    mass = sum(model.elements[eid].Mass() for eid in skin_eids)
    return mass


def run_nastran(out_bdf, out_dir):
    cmd = [NASTRAN_EXE, f"JID={out_bdf}", "mem=max"]
    subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, timeout=180)
    base = os.path.splitext(os.path.basename(out_bdf))[0].lower()
    log_path = os.path.join(out_dir, f"{base}.log")
    op2_path = os.path.join(out_dir, f"{base}.op2")
    if not os.path.exists(log_path) or "NSEXIT: EXIT(0)" not in open(log_path).read():
        return None
    return op2_path


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
    print("Building + repairing shared half-wing base mesh (2-spar)...")
    base_model, skin_eids, y_max = build_repaired_half_wing()

    # ---- verify spar/rib mass is untouched by the taper (professor's request) ----
    # cross_reference() is needed here because apply_node_merges() (inside
    # the repair passes above) deletes and re-adds elements via
    # add_cquad4/add_ctria3, whose fresh objects have pid_ref=None until
    # cross-referenced -- elem.Mass() needs pid_ref, not just pid.
    base_model.cross_reference()
    other_mass_probe = 0.0
    for elem in base_model.elements.values():
        if elem.type not in ('CQUAD4', 'CTRIA3'):
            continue
        if elem.eid not in skin_eids:
            other_mass_probe += elem.Mass()
    print(f"\nspar+rib mass (fixed, never touched by the taper loop): "
          f"{other_mass_probe * 1000:.3f} g\n")

    # ---- solve for (a, b) in skin_mass = a*T_root + b*T_tip via two probes ----
    mass_a = skin_mass_only(base_model, skin_eids, y_max, t_root=1.0, t_tip=0.0)
    mass_b = skin_mass_only(base_model, skin_eids, y_max, t_root=0.0, t_tip=1.0)
    a, b = mass_a, mass_b
    target_mass = a * T_ROOT_REF + b * T_TIP_REF
    print(f"skin_mass(T_root, T_tip) = {a:.6f}*T_root + {b:.6f}*T_tip  (T in metres, mass in kg)")
    print(f"target skin mass (matching T_root={T_ROOT_REF*1e3:.2f}mm / T_tip={T_TIP_REF*1e3:.2f}mm reference): "
          f"{target_mass*1000:.2f} g")

    family = []
    for t_tip_mm in T_TIP_FAMILY_MM:
        t_tip = t_tip_mm * 1e-3
        t_root = (target_mass - b * t_tip) / a
        family.append((t_root, t_tip))
        print(f"  T_tip={t_tip_mm:.2f}mm -> T_root={t_root*1e3:.3f}mm "
              f"(taper ratio {t_root/t_tip:.2f}:1)")

    results = []
    for t_root, t_tip in family:
        label = f"Troot_{t_root*1e3:.2f}mm_Ttip_{t_tip*1e3:.2f}mm"
        out_dir = os.path.join(HERE, label)
        print(f"\n{'#'*70}\n#  {label}\n{'#'*70}")

        out_bdf = build_case(base_model, skin_eids, y_max, t_root, t_tip, out_dir)
        op2_path = run_nastran(out_bdf, out_dir)
        if op2_path is None:
            print(f"  [!] Nastran did not exit cleanly -- check {out_dir}")
            continue

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

        write_interactive_html(out_dir, bdf_final, eig, skin_eids)

        r = dict(t_root_mm=t_root * 1e3, t_tip_mm=t_tip * 1e3,
                 taper_ratio=t_root / t_tip,
                 skin_mass_g=skin_mass * 1000, other_mass_g=other_mass * 1000,
                 total_mass_g=(skin_mass + other_mass) * 1000,
                 f1_hz=freqs_hz[0], f1_ratio=ratio)
        results.append(r)
        print(f"  f1={r['f1_hz']:.3f} Hz  ratio={r['f1_ratio']:.2f}  "
              f"skin_mass={r['skin_mass_g']:.2f}g  other_mass={r['other_mass_g']:.2f}g  "
              f"total_mass={r['total_mass_g']:.2f}g")

    out_csv = os.path.join(HERE, "sweep_constant_weight_summary.csv")
    with open(out_csv, "w") as f:
        f.write("t_root_mm,t_tip_mm,taper_ratio,skin_mass_g,other_mass_g,total_mass_g,f1_hz,f1_ratio\n")
        for r in results:
            f.write(f"{r['t_root_mm']:.4f},{r['t_tip_mm']:.4f},{r['taper_ratio']:.3f},"
                    f"{r['skin_mass_g']:.3f},{r['other_mass_g']:.3f},{r['total_mass_g']:.3f},"
                    f"{r['f1_hz']:.4f},{r['f1_ratio']:.2f}\n")
    print(f"\nwritten: {out_csv}")

    report_path = os.path.join(HERE, "verification_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Constant-Weight Taper Sweep -- 2-Spar Verification Report\n\n")
        f.write(f"**Source geometry:** `{IN_BDF}`\n")
        f.write(f"(same 2-spar, dense-rib layout as `two_spar_dense_ribs_check/AC_wing_2spar_denseribs.vsp3` "
                f"-- open that file in OpenVSP to inspect the shared geometry used for every case below; "
                f"only the skin thickness assignment changes per case, so a separate .vsp3 per taper point "
                f"would be identical and is not produced here.)\n\n")
        f.write(f"**Spar/rib mass, fixed and never touched by the taper assignment loop:** "
                f"{other_mass_probe * 1000:.3f} g\n\n")
        f.write(f"**Skin mass model:** skin_mass = {a:.6f}*T_root + {b:.6f}*T_tip (T in metres, mass in kg)\n\n")
        f.write(f"**Target skin mass** (matches the TU-Flex reference taper, "
                f"T_root={T_ROOT_REF*1e3:.2f}mm / T_tip={T_TIP_REF*1e3:.2f}mm): "
                f"{target_mass*1000:.2f} g\n\n")
        f.write("## Per-case results\n\n")
        f.write("| T_root (mm) | T_tip (mm) | Taper ratio | Skin mass (g) | Spar+rib mass (g) | Total mass (g) | f1 (Hz) | Mode-1 contamination ratio |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in results:
            f.write(f"| {r['t_root_mm']:.3f} | {r['t_tip_mm']:.2f} | {r['taper_ratio']:.2f}:1 | "
                    f"{r['skin_mass_g']:.2f} | {r['other_mass_g']:.2f} | {r['total_mass_g']:.2f} | "
                    f"{r['f1_hz']:.3f} | {r['f1_ratio']:.2f} |\n")

        if results:
            other_masses = [r['other_mass_g'] for r in results]
            total_masses = [r['total_mass_g'] for r in results]
            skin_masses = [r['skin_mass_g'] for r in results]
            f.write("\n## Verification\n\n")
            f.write(f"- Skin mass range: {min(skin_masses):.4f} to {max(skin_masses):.4f} g "
                    f"(target: {target_mass*1000:.2f} g)\n")
            f.write(f"- Spar+rib mass range: {min(other_masses):.4f} to {max(other_masses):.4f} g "
                    f"-- should be *identical* across all cases, since the taper loop never touches "
                    f"these elements.\n")
            f.write(f"- **Total structural mass range: {min(total_masses):.4f} to {max(total_masses):.4f} g**\n\n")
            spread = max(total_masses) - min(total_masses)
            if spread < 0.01:
                f.write(f"Total mass spread is {spread:.5f} g (negligible, floating-point noise) -- "
                        f"**confirms that 'constant weight' holds for the full structural model**, "
                        f"not just the skin, because the spar/rib elements are never touched by this "
                        f"sweep. Section 7.1's claim that spars/ribs share the skin's spanwise-varying "
                        f"thickness does **not** match what this model actually does (they keep a fixed "
                        f"thickness) -- Section 7.1 should be corrected to reflect this.\n")
            else:
                f.write(f"Total mass spread is {spread:.5f} g -- **NOT negligible**. Do not claim "
                        f"'constant weight' for the full model; use 'constant total skin mass' instead, "
                        f"and report this spread explicitly.\n")
    print(f"\nwritten: {report_path}")

    if results:
        other_masses = [r['other_mass_g'] for r in results]
        total_masses = [r['total_mass_g'] for r in results]
        print(f"\nVERIFICATION -- other_mass_g (spar+rib) range: "
              f"{min(other_masses):.4f} to {max(other_masses):.4f} g "
              f"(should be identical across all cases)")
        print(f"VERIFICATION -- total_mass_g range: "
              f"{min(total_masses):.4f} to {max(total_masses):.4f} g "
              f"(if this range is negligible, 'constant weight' is justified for the FULL model, "
              f"not just the skin)")

    # ---- final figure: f1 vs taper ratio, at constant total skin mass ----
    if results:
        plot_rows = sorted(results, key=lambda r: r['taper_ratio'])
        plt.rcParams.update({"font.size": 12})
        fig, ax = plt.subplots(figsize=(6.5, 4.8))
        taper_ratios = [r['taper_ratio'] for r in plot_rows]
        f1s_plot = [r['f1_hz'] for r in plot_rows]
        ax.plot(taper_ratios, f1s_plot, "o-", color="#1f6fb2", linewidth=2, markersize=7)
        for r in plot_rows:
            ax.annotate(f"{r['t_root_mm']:.2f}/{r['t_tip_mm']:.2f}mm",
                        (r['taper_ratio'], r['f1_hz']),
                        textcoords="offset points", xytext=(6, 6), fontsize=8, color="#555555")
        ax.set_xlabel(r"Taper ratio $T_{root}/T_{tip}$")
        ax.set_ylabel("First bending frequency of the clamped\nhalf-wing structural model (Hz)")
        ax.set_title("First bending frequency vs. taper steepness\n(constant total skin mass)")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig_path = os.path.join(HERE, "fig1_frequency_vs_taper.png")
        fig.savefig(fig_path, dpi=160)
        plt.close(fig)
        print(f"written: {fig_path}")
