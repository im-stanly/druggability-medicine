"""Hypothesis: Simple pocket descriptors (volume, hydrophobicity, polarity)
classify druggable vs non-druggable pockets with AUROC > 0.75,
on both PDB and AlphaFold structures.

Results summary (4138 pockets, 438 proteins):
  - 8 simple features: atom count, residue count, density,
    hydrophobic/polar/charged/positive/negative fractions
  - Label: druggable if pocket center ≤ 6Å from any non-solvent ligand
  - RandomForest (100 trees, max_depth=5)
  - Stratified 5-fold CV by protein ID

Usage:
    uv run python verify_simple_classifier.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "src"))

# ── config ────────────────────────────────────────────────────────────
POCKETS_CSV = Path("data/08_reporting/p2rank_pockets.csv")
LIGAND_RADIUS = 6.0
FEATURE_RADIUS = 6.0

HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO", "PHE", "TRP", "CYS", "GLY"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "HIS"}
POSITIVE = {"LYS", "ARG"}
NEGATIVE = {"ASP", "GLU"}

SKIP_LIGANDS = {"HOH", "DOD", "WAT", "GOL", "EDO", "SO4", "PO4", "ACT", "PEG",
                "BME", "DMS", "FMT", "EPE", "CIT", "TRS", "MPD", "MG", "CL",
                "NA", "K", "CA", "ZN", "MN", "FE", "CO", "NI", "CD", "HG"}


def main():
    import gzip, io
    from Bio.PDB import MMCIFParser
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score, StratifiedGroupKFold
    from sklearn.metrics import roc_auc_score

    # ── 1. Featurize all pockets ────────────────────────────────────
    print("Featurizing pockets...")
    df = pd.read_csv(POCKETS_CSV)
    structures = {}

    features_list = []
    labels = []
    pdb_ids = []

    for i, (_, row) in enumerate(df.iterrows()):
        pdb_id = row["pdb_id"]
        if i % 500 == 0:
            print(f"  {i}/{len(df)}")

        if pdb_id not in structures:
            cif_path = Path(f"data/scrapped/structures/{pdb_id}/PDB-{pdb_id}.cif.gz")
            if not cif_path.exists():
                continue
            with gzip.open(cif_path, "rt", encoding="utf-8") as f:
                content = f.read()
            parser = MMCIFParser(QUIET=True)
            try:
                structures[pdb_id] = parser.get_structure("s", io.StringIO(content))
            except Exception:
                continue

        s = structures[pdb_id]
        center = np.array([float(row["center_x"]), float(row["center_y"]), float(row["center_z"])])

        # Featurize
        atoms = list(s.get_atoms())
        coords = np.array([a.coord for a in atoms])
        dists = np.linalg.norm(coords - center[None, :], axis=1)
        mask = dists <= FEATURE_RADIUS
        if mask.sum() < 5:
            continue

        near_atoms = [a for a, m in zip(atoms, mask) if m]
        seen = set()
        h = p = pos = neg = 0
        for a in near_atoms:
            res = a.get_parent()
            if res.id[0].strip():
                continue
            rn = res.resname.strip().upper()
            rid = res.id
            key = (rn, rid)
            if key in seen:
                continue
            seen.add(key)
            if rn in HYDROPHOBIC: h += 1
            elif rn in POLAR: p += 1
            elif rn in POSITIVE: pos += 1
            elif rn in NEGATIVE: neg += 1

        n_atoms = int(mask.sum())
        n_res = max(len(seen), 1)

        features_list.append({
            "n_atoms": n_atoms,
            "n_residues": n_res,
            "atom_density": n_atoms / (4/3 * np.pi * FEATURE_RADIUS**3),
            "frac_hydrophobic": h / n_res,
            "frac_polar": p / n_res,
            "frac_charged": (pos + neg) / n_res,
            "frac_positive": pos / n_res,
            "frac_negative": neg / n_res,
        })

        # Label
        label = 0
        for res in s.get_residues():
            if not res.id[0].strip():
                continue
            if res.resname.strip() in SKIP_LIGANDS:
                continue
            for atom in res:
                if atom.element == "H":
                    continue
                if np.linalg.norm(np.array(atom.coord) - center) <= LIGAND_RADIUS:
                    label = 1
                    break
            if label == 1:
                break

        labels.append(label)
        pdb_ids.append(pdb_id)

    X = pd.DataFrame.from_records(features_list)
    y = pd.Series(labels, name="druggable")
    groups = pd.Series(pdb_ids, name="pdb_id")

    print(f"\n  Total: {len(X)} pockets, {y.sum()} druggable ({100*y.mean():.1f}%)")
    print(f"  Proteins: {groups.nunique()}")

    # ── 2. Cross-validate ──────────────────────────────────────────
    print(f"\n{'='*60}")
    print("5-fold CV (stratified by protein):")
    print(f"  Features: {list(X.columns)}")

    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    model = RandomForestClassifier(
        n_estimators=100, max_depth=5, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )

    aurocs = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups=groups)):
        X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]

        model.fit(X_tr, y_tr)
        y_prob = model.predict_proba(X_te)[:, 1]
        auc = roc_auc_score(y_te, y_prob)
        aurocs.append(auc)
        print(f"  Fold {fold+1}: AUROC = {auc:.3f}  (n={len(X_te)}, druggable={y_te.sum()})")

    mean_auc = np.mean(aurocs)
    std_auc = np.std(aurocs)
    print(f"\n  Mean AUROC: {mean_auc:.3f} ± {std_auc:.3f}")
    print(f"  {'✓ > 0.75' if mean_auc > 0.75 else '✗ < 0.75'}")

    # ── 3. Final model + feature importances ────────────────────────
    model.fit(X, y)
    print(f"\n  Feature importances:")
    for name, imp in sorted(zip(X.columns, model.feature_importances_), key=lambda x: -x[1]):
        print(f"    {name:<20s} {imp:.4f}")

    # ── 4. Save ────────────────────────────────────────────────────
    result = {
        "hypothesis": "Simple pocket descriptors classify druggable vs non-druggable with AUROC > 0.75",
        "auroc_mean": float(mean_auc),
        "auroc_std": float(std_auc),
        "auroc_folds": [float(a) for a in aurocs],
        "passed": bool(mean_auc > 0.75),
        "n_pockets": int(len(X)),
        "n_druggable": int(y.sum()),
        "n_proteins": int(groups.nunique()),
        "features": list(X.columns),
        "feature_importances": {str(name): float(c) for name, c in zip(X.columns, model.feature_importances_)},
        "model": "RandomForest(n_estimators=100, max_depth=5, class_weight=balanced)",
        "labeling": f"druggable if pocket center ≤ {LIGAND_RADIUS}Å from non-solvent ligand atom",
    }
    out_path = Path("data/08_reporting/simple_classifier_results.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
