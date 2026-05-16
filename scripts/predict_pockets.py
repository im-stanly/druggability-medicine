#!/usr/bin/env python3
"""Predict binding pockets on a protein structure.

Usage:
    uv run python scripts/predict_pockets.py --cif-path data/.../PDB-XXXX.cif.gz
    uv run python scripts/predict_pockets.py --all-af
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import KDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from druggability.pocket_mining.parser import parse_cif
from druggability.pocket_mining.surface import generate_surface_points
from druggability.pocket_mining.features import featurize
from druggability.pocket_mining.model import PocketClassifier
from druggability.pocket_mining.constants import (
    PROBABILITY_THRESHOLD, CLUSTER_RADIUS, MIN_CLUSTER_SIZE,
)


def find_pockets(points, probs, threshold=PROBABILITY_THRESHOLD,
                 radius=CLUSTER_RADIUS, min_size=MIN_CLUSTER_SIZE):
    """Group high-probability points into pocket regions."""
    mask = probs >= threshold
    if mask.sum() < min_size:
        return []

    pts = points[mask]
    pr = probs[mask]
    order = np.argsort(pr)[::-1]
    pts, pr = pts[order], pr[order]

    tree = KDTree(pts)
    taken = np.zeros(len(pts), dtype=bool)
    pockets = []

    for i in range(len(pts)):
        if taken[i]:
            continue
        nb = [j for j in tree.query_ball_point(pts[i], radius) if not taken[j]]
        if len(nb) < min_size:
            continue
        for j in nb:
            taken[j] = True

        pocket_pts = pts[nb]
        pockets.append({
            "rank": len(pockets) + 1,
            "center": pocket_pts.mean(axis=0),
            "score": float(pr[nb].mean()),
            "n": len(nb),
            "points": pocket_pts,
        })

    return pockets


def predict_one(cif_path, model_path="models/pocket_classifier.pkl",
                n_pts=3000, threshold=PROBABILITY_THRESHOLD,
                out_dir=None):
    """Run prediction on a single .cif.gz file."""
    cif_path = Path(cif_path)
    clf = PocketClassifier.load(model_path)
    parsed = parse_cif(cif_path)

    print(f"Protein: {len(parsed.protein.coords)} atoms, "
          f"{len(parsed.ligands)} ligands (ground truth)")

    pts = generate_surface_points(parsed.protein, n_points=n_pts)
    X = featurize(pts, parsed.protein)
    probs = clf.predict_proba(X)

    n_hit = (probs >= threshold).sum()
    print(f"Surface points: {len(pts)}, "
          f"{n_hit} above threshold ({100*n_hit/len(pts):.1f}%)")

    pockets = find_pockets(pts, probs, threshold)
    print(f"Pocket regions found: {len(pockets)}")

    # --- save ---
    out_dir = Path(out_dir or cif_path.parent)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = cif_path.stem.replace(".cif", "")

    # summary CSV
    with open(out_dir / f"{base}_pockets.csv", "w") as f:
        f.write("rank,score,cx,cy,cz,n_points\n")
        for p in pockets:
            cx, cy, cz = p["center"]
            f.write(f"{p['rank']},{p['score']:.4f},{cx:.3f},{cy:.3f},{cz:.3f},{p['n']}\n")

    # point-level CSV
    with open(out_dir / f"{base}_points.csv", "w") as f:
        f.write("x,y,z,probability\n")
        for pt, prob in zip(pts, probs):
            f.write(f"{pt[0]:.3f},{pt[1]:.3f},{pt[2]:.3f},{prob:.4f}\n")

    # pocket PDBs for PyMOL
    for p in pockets:
        pdb = out_dir / f"{base}_pocket{p['rank']}.pdb"
        with open(pdb, "w") as f:
            for j, pt in enumerate(p["points"]):
                f.write(f"HETATM{j+1:5d}  C   POC P{p['rank']:1d}{j+1:4d}    "
                        f"{pt[0]:8.3f}{pt[1]:8.3f}{pt[2]:8.3f}  1.00  0.00      C\n")
            f.write("END\n")

    print(f"Output → {out_dir}/")
    print(f"  {base}_pockets.csv")
    print(f"  {base}_points.csv")
    for p in pockets:
        print(f"  {base}_pocket{p['rank']}.pdb")

    # compare with ground truth
    if parsed.ligands:
        print(f"\nGround truth ligands:")
        for lig in parsed.ligands:
            print(f"  {lig.name:5s} at ({lig.center[0]:.1f}, {lig.center[1]:.1f}, {lig.center[2]:.1f})")
        for p in pockets:
            for lig in parsed.ligands:
                d = np.linalg.norm(p["center"] - lig.center)
                if d < 8.0:
                    print(f"  pocket {p['rank']} → {d:.1f} Å from {lig.name}")
                    break

    return pockets


def main():
    ap = argparse.ArgumentParser(description="Predict binding pockets")
    ap.add_argument("--cif-path", "-p", help="Path to .cif.gz file")
    ap.add_argument("--all-af", action="store_true",
                    help="Run on all AlphaFold structures in data/")
    ap.add_argument("--model-path", "-m", default="models/pocket_classifier.pkl")
    ap.add_argument("--n-points", type=int, default=3000)
    ap.add_argument("--threshold", type=float, default=PROBABILITY_THRESHOLD)
    args = ap.parse_args()

    if args.all_af:
        data = Path("data/01_raw/proteins")
        for af in sorted(data.glob("*/AF-*.cif.gz")):
            print(f"\n{'='*50}\n{af.parent.name}/{af.name}\n{'='*50}")
            try:
                predict_one(af, args.model_path, args.n_points, args.threshold)
            except Exception as e:
                print(f"ERROR: {e}")
    elif args.cif_path:
        predict_one(args.cif_path, args.model_path, args.n_points, args.threshold)
    else:
        print("Pass --cif-path or --all-af")
        sys.exit(1)


if __name__ == "__main__":
    main()
