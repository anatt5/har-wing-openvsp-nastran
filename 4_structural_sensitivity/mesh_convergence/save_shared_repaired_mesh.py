"""
Save the repaired-but-not-yet-taper/material-assigned half-wing mesh for
BaseLen=0.50 as a standalone canonical snapshot, so other studies (the
taper sweep, the AR sweep) can read this EXACT repaired mesh directly
instead of running their own separate repair pipeline on the raw
geometry. This removes any possibility of two independently-implemented
repair pipelines producing subtly different meshes from the same raw
input -- which is exactly what caused a real discrepancy between the
taper-sweep and mesh-convergence results earlier in this thesis (traced
back to a pyNastran `pid_ref` staleness bug, not to mesh non-determinism).

Depends on run_mesh_convergence.py having already been run once for
BaseLen=0.50 (reads its output "BaseLen_0.50/AC_wing.bdf", and reuses its
mesh-repair functions directly via import).

Run:  python save_shared_repaired_mesh.py
"""
import os
from pyNastran.bdf.bdf import BDF
import run_mesh_convergence as mc

HERE = os.path.dirname(os.path.abspath(__file__))
raw_bdf = os.path.join(HERE, "BaseLen_0.50", "AC_wing.bdf")
raw_bdf_bdf = raw_bdf + ".bdf"

model = BDF(debug=False)
model.read_bdf(raw_bdf_bdf, xref=True, validate=False)

# ---- same repair pipeline used by every other structural script in this
# thesis (see run_mesh_convergence.py for the canonical implementation):
# fix pre-existing duplicate-ID quads from the raw OpenVSP mesh, then
# equivalence near-duplicate nodes, then iterate degenerate-quad repair
# and weak-connectivity-node merging until both converge.
mc.fix_preexisting_duplicate_id_quads(model)
mc.equivalence_duplicate_nodes(model, tol=mc.NODE_EQUIV_TOL_MM * 1e-3)
for _ in range(6):
    n_angle = mc.repair_degenerate_quads(model)
    n_weak = mc.merge_weak_connectivity_nodes(model)
    if n_angle == 0 and n_weak == 0:
        break

# ---- keep only the half-wing (y <= 0): OpenVSP mirrors the wing about
# the root plane, and the structural studies only ever clamp/analyze one
# half.
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

# ---- any node that still has too little connectivity after the repair
# loop above is rigidized to its best neighbor with an RBE2, rather than
# left as a numerically weak point in the eigenvalue solve.
mc.rigidize_remaining_weak_nodes(model)
model.cross_reference()

out_path = os.path.join(HERE, "BaseLen_0.50", "AC_wing_repaired_half_wing.bdf")
model.write_bdf(out_path)
print(f"written: {out_path}")
print(f"nodes={len(model.nodes)} elements={len(model.elements)}")
