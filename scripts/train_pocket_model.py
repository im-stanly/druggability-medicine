#!/usr/bin/env python3
"""Train a pocket classifier on PDB structures with bound ligands.

Scans data/01_raw/proteins/ for PDB .cif.gz files, extracts protein +
ligand info, generates surface points, labels them by distance from
the ligand, featurises, and trains a Random Forest.

Usage:
    uv run python scripts/train_pocket_model.py
"""

import sys
import time
from pathlib import Path

import numpy as np
from sklearn.model_selection import cross_val_score, StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from druggability.pocket_mining.parser import parse_cif, get_ligand_atoms
from druggability.pocket_mining.surface import generate_surface_points, label_points
from druggability.pocket_mining.features import featurize, feature_names
from druggability.pocket_mining.model import PocketClassifier
from druggability.pocket_mining.constants import POCKET_RADIUS, NON_POCKET_RADIUS


def collect_data(data_dir: Path, n_pts: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """Walk protein dirs, extract features and labels from each PDB structure."""
    X_parts, y_parts = [], []
    ok, skip = 0, 0

    dirs = sorted(p for p in data_dir.iterdir() if p.is_dir() and not p.name.startswith("."))
    print(f"Found {len(dirs)} protein directories\n")

    for i, d in enumerate(dirs):
        pdb = list(d.glob("PDB-*.cif.gz"))
        if not pdb:
            skip += 1; continue

        try:
            t0 = time.time()
            parsed = parse_cif(pdb[0])
            pts = generate_surface_points(parsed.protein, n_points=n_pts)
            lab = label_points(pts, get_ligand_atoms(parsed),
                               POCKET_RADIUS, NON_POCKET_RADIUS)
            valid = lab >= 0
            if valid.sum() == 0:
                skip += 1; continue

            X_parts.append(featurize(pts[valid], parsed.protein))
            y_parts.append(lab[valid])
            ok += 1

            n_pos = (lab[valid] == 1).sum()
            print(f"  [{i+1:2d}/{len(dirs)}] {parsed.pdb_id:6s}  "
                  f"{valid.sum():5d} pts  +{n_pos:4d}/-{valid.sum()-n_pos:4d}  "
                  f"{len(parsed.ligands)} ligands  {time.time()-t0:.1f}s")
        except Exception as e:
            skip += 1
            print(f"  [{i+1:2d}/{len(dirs)}] {d.name}: {e}")

    X = np.vstack(X_parts)
    y = np.concatenate(y_parts)
    print(f"\n{ok} proteins, {skip} skipped, {len(X)} points "
          f"(+{(y==1).sum()}/-{(y==0).sum()})\n")
    return X, y


def main():
    print("=" * 55)
    print("POCKET CLASSIFIER  (15 features, Random Forest)")
    print("=" * 55)

    # 1. Collect
    print("\n[1/3] Collecting training data ...")
    X, y = collect_data(Path("data/01_raw/proteins"))

    # 2. Train
    print("[2/3] Training ...")
    t0 = time.time()
    clf = PocketClassifier()
    clf.fit(X, y, feature_names=feature_names())
    print(f"      done in {time.time()-t0:.1f}s")

    scores = clf.score(X, y)
    print(f"      train acc={scores['accuracy']:.3f}  "
          f"roc={scores['roc_auc']:.3f}  f1={scores['f1']:.3f}")

    # 3. Cross-validate
    print("[3/3] 5-fold CV ...")
    cv = StratifiedKFold(5, shuffle=True, random_state=42)
    for metric in ("accuracy", "roc_auc", "f1"):
        s = cross_val_score(clf.model, clf.scaler.transform(X), y,
                            cv=cv, scoring=metric)
        print(f"      cv {metric:10s}: {s.mean():.3f} ± {s.std():.3f}")

    # Save
    path = Path("models/pocket_classifier.pkl")
    clf.save(path)
    print(f"\nSaved → {path.resolve()}")

    # Top features
    print("\nFeature importances:")
    for f in clf.importances()[:10]:
        print(f"  {f['feature']:20s} {f['importance']:.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
