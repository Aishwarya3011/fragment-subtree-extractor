#!/usr/bin/env python3
import os
import argparse
import csv
import math
from typing import Dict, List, Tuple, Optional
from collections import defaultdict, deque

def parse_giant_swc(path: str):
    nodes: Dict[int, Tuple[float, float, float, float, int, int]] = {}
    children: Dict[int, List[int]] = defaultdict(list)
    roots: List[int] = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            nid = int(parts[0])
            ntype = int(float(parts[1]))
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            r = float(parts[5])
            parent = int(parts[6])

            nodes[nid] = (x, y, z, r, ntype, parent)

    for nid, (x, y, z, r, ntype, parent) in nodes.items():
        if parent == -1:
            # Start of a fragement. Potential root.
            roots.append(nid)
        else:
            if parent in nodes:
                children[parent].append(nid)
            else:
                # No parent. Could still be root.
                roots.append(nid)

    roots = sorted(set(roots))
    return nodes, children, roots

def extract_subtree(children: Dict[int, List[int]], root: int) -> List[int]:
    # BFS from every child
    out = []
    q = deque([root])
    seen = set()
    while q:
        u = q.popleft()
        if u in seen:
            continue
        seen.add(u)
        out.append(u)
        for v in children.get(u, []):
            if v not in seen:
                q.append(v)
    return out

def write_swc_for_root(
    out_path: str,
    nodes,
    children,
    root_id: int,
    root_type: int = 1,
    neurite_type: int = 0,
):
    subtree = extract_subtree(children, root_id)

    order = []
    q = deque([root_id])
    seen = set()
    while q:
        u = q.popleft()
        if u in seen:
            continue
        seen.add(u)
        order.append(u)
        for v in sorted(children.get(u, [])):
            if v not in seen:
                q.append(v)

    for nid in subtree:
        if nid not in seen:
            order.append(nid)

    new_id = {old: i + 1 for i, old in enumerate(order)}

    with open(out_path, "w") as f:
        f.write("# id type x y z radius parent\n")

        for old in order:
            x, y, z, r, _, parent = nodes[old]
            swc_type = root_type if old == root_id else neurite_type

            if old == root_id:
                swc_parent = -1
            else:
                if parent not in new_id:
                    swc_parent = new_id[root_id]
                else:
                    swc_parent = new_id[parent]

            f.write(
                f"{new_id[old]} {swc_type} {x:.4f} {y:.4f} {z:.4f} {r:.4f} {swc_parent}\n"
            )

def read_anchors_csv(path: str) -> List[Tuple[float, float, float]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("anchors.csv has no header row.")

        fields = list(reader.fieldnames)
        lower = [c.strip().lower() for c in fields]

        def find_col(target: str) -> Optional[str]:
            for orig, lo in zip(fields, lower):
                if lo == target:
                    return orig
            for orig, lo in zip(fields, lower):
                if lo.endswith(f"_{target}") or lo.endswith(target) or f"_{target}_" in lo:
                    if target in lo:
                        return orig
            for orig, lo in zip(fields, lower):
                if target in lo:
                    return orig
            return None

        xcol = find_col("x")
        ycol = find_col("y")
        zcol = find_col("z")

        anchors: List[Tuple[float, float, float]] = []
        for row in reader:
            x = float(row[xcol])
            y = float(row[ycol])
            z = float(row[zcol])
            anchors.append((x, y, z))

    return anchors

def in_cube(px: float, py: float, pz: float, cx: float, cy: float, cz: float, half: float) -> bool:
    return (abs(px - cx) <= half) and (abs(py - cy) <= half) and (abs(pz - cz) <= half)

def anchor_swc_filename(sx: float, sy: float, sz: float) -> str:
    sx_i = int(round(sx))
    sy_i = int(round(sy))
    sz_i = int(round(sz))
    return f"{sx_i}x-{sy_i}y-{sz_i}z.swc"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--giant_swc", required=True, help="Path to raw reconstructions (giant SWC)")
    ap.add_argument("--anchors_csv", required=True, help="CSV with anchor coordinates (manually annotated)")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--min_nodes", type=int, default=10, help="Skip roots with fewer than this many nodes")
    ap.add_argument("--cube_half", type=float, required=True, help="Half-width of the cube around each anchor")

    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    anchors = read_anchors_csv(args.anchors_csv)
    if not anchors:
        raise ValueError("No anchors found")
    print(f"Loading {len(anchors)} anchors")

    anchor_dirs: List[str] = []
    for i, (sx, sy, sz) in enumerate(anchors):
        name = anchor_swc_filename(sx, sy, sz)
        d = os.path.join(args.out_dir, name)
        os.makedirs(d, exist_ok=True)
        anchor_dirs.append(d)

    nodes, children, roots = parse_giant_swc(args.giant_swc)

    # For printing progress
    assigned = 0
    for root in roots:
        subtree = extract_subtree(children, root)
        # To curb stray fraegments. May remove later.
        if len(subtree) < args.min_nodes:
            continue

        rx, ry, rz, *_ = nodes[root]
        candidate_anchors: List[int] = []

        for si, (sx, sy, sz) in enumerate(anchors):
            hit = False
            for nid in subtree:
                x, y, z, *_ = nodes[nid]
                if in_cube(x, y, z, sx, sy, sz, args.cube_half):
                    hit = True
                    break
            if hit:
                candidate_anchors.append(si)

        if not candidate_anchors:
            continue

        for si in candidate_anchors:
            out_path = os.path.join(
            anchor_dirs[si],
            f"root_{root:08d}.swc",
        )
            write_swc_for_root(out_path, nodes, children, root)
            assigned += 1

        if assigned % 100 == 0:
            print(f"Written roots: {assigned}")

if __name__ == "__main__":
    main()
