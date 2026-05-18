"""Kedro nodes for the pocket model pipeline.

Each node is a pure function: takes inputs, returns outputs.
Kedro wires them together based on the catalog and pipeline definition.
"""

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from druggability.pocket_mining.parser import parse_cif, get_ligand_atoms
from druggability.pocket_mining.surface import generate_surface_points, label_points
from druggability.pocket_mining.features import featurize, feature_names
from druggability.pocket_mining.model import PocketClassifier
from druggability.pocket_mining.constants import POCKET_RADIUS, NON_POCKET_RADIUS

logger = logging.getLogger(__name__)


def train_pocket_classifier(data_dir: str, model_path: str,
                            n_points: int = 2000) -> dict:
    """Collect training data from protein directories, train, save, return metrics.

    Args:
        data_dir: Path to directory with protein subdirectories (each holding
                  PDB-*.cif.gz files).
        model_path: Where to save the trained .pkl model.
        n_points: Surface points to generate per protein.

    Returns:
        Dict with training metrics (accuracy, roc_auc, f1, cv scores, n_proteins, n_points).
    """
    data_dir = Path(data_dir)
    logger.info("Collecting training data from %s ...", data_dir)

    X_parts, y_parts = [], []
    n_proteins, n_skipped = 0, 0

    for d in sorted(p for p in data_dir.iterdir()
                    if p.is_dir() and not p.name.startswith(".")):
        pdb = list(d.glob("PDB-*.cif.gz"))
        if not pdb:
            n_skipped += 1
            continue
        try:
            t0 = time.time()
            parsed = parse_cif(pdb[0])
            pts = generate_surface_points(parsed.protein, n_points=n_points)
            lab = label_points(pts, get_ligand_atoms(parsed),
                               POCKET_RADIUS, NON_POCKET_RADIUS)
            valid = lab >= 0
            if valid.sum() == 0:
                n_skipped += 1
                continue
            X_parts.append(featurize(pts[valid], parsed.protein))
            y_parts.append(lab[valid])
            n_proteins += 1
            n_pos = (lab[valid] == 1).sum()
            logger.debug("  %s: %d pts (+%d/-%d) %.1fs",
                         parsed.pdb_id, valid.sum(), n_pos,
                         valid.sum() - n_pos, time.time() - t0)
        except Exception as e:
            n_skipped += 1
            logger.warning("  Skipping %s: %s", d.name, e)

    if not X_parts:
        raise RuntimeError("No training data collected from %s" % data_dir)

    X = np.vstack(X_parts)
    y = np.concatenate(y_parts)

    logger.info("Collected %d proteins (%d skipped), %d points (+%d/-%d)",
                n_proteins, n_skipped, len(X), (y == 1).sum(), (y == 0).sum())

    # ── train ──────────────────────────────────────────────────────
    logger.info("Training Random Forest ...")
    t0 = time.time()
    clf = PocketClassifier()
    clf.fit(X, y, feature_names=feature_names())
    train_time = time.time() - t0

    train_scores = clf.score(X, y)
    logger.info("Train: acc=%.3f  roc=%.3f  f1=%.3f  (%.1fs)",
                train_scores["accuracy"], train_scores["roc_auc"],
                train_scores["f1"], train_time)

    # ── cross-validate ─────────────────────────────────────────────
    logger.info("5-fold CV ...")
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    cv = StratifiedKFold(5, shuffle=True, random_state=42)

    cv_results = {}
    for metric in ("accuracy", "roc_auc", "f1"):
        scores = cross_val_score(clf.model, Xs, y, cv=cv, scoring=metric)
        cv_results[f"cv_{metric}_mean"] = float(scores.mean())
        cv_results[f"cv_{metric}_std"] = float(scores.std())
        logger.info("  cv %s: %.3f ± %.3f", metric, scores.mean(), scores.std())

    # ── save model ─────────────────────────────────────────────────
    clf.save(model_path)
    logger.info("Model saved → %s", model_path)

    # ── top features ───────────────────────────────────────────────
    top = clf.importances()[:5]
    logger.info("Top features: %s",
                ", ".join(f"{f['feature']}({f['importance']:.3f})" for f in top))

    return {
        "n_proteins": n_proteins,
        "n_skipped": n_skipped,
        "n_points": len(X),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "pos_ratio": float((y == 1).sum() / len(y)),
        "train_accuracy": train_scores["accuracy"],
        "train_roc_auc": train_scores["roc_auc"],
        "train_f1": train_scores["f1"],
        "train_time_s": train_time,
        **cv_results,
    }
