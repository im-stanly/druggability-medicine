from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PocketLabelingParams:
    max_ligand_distance_angstrom: float = 6.0


def label_pockets_by_ligand_distance(
    pockets: pd.DataFrame,
    protein_ligand_ds: Iterable[dict[str, Any]],
    params: dict[str, Any],
) -> pd.DataFrame:
    """Assign a weak label to each pocket based on proximity to any ligand centroid.

    Labeling rule (baseline):
      * druggable = 1 if min distance from pocket center to any ligand centroid <= threshold
      * druggable = 0 otherwise

    This is a baseline heuristic to get a trainable dataset.
    """
    p = PocketLabelingParams(**params)

    # Build pdb_id -> ligand centroids map from dataset items.
    lig_map: dict[str, list[tuple[float, float, float]]] = {}
    for item in protein_ligand_ds:
        pdb_id = item["pdb_id"]
        centroids = [tuple(v["centroid"]) for v in item.get("ligands", {}).values()]
        lig_map[pdb_id] = centroids

    df = pockets.copy()
    for c in ["center_x", "center_y", "center_z"]:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in pockets dataframe")

    min_dists: list[float | None] = []
    druggable: list[int] = []

    for _, row in df.iterrows():
        pdb_id = row["pdb_id"]
        cx, cy, cz = float(row["center_x"]), float(row["center_y"]), float(row["center_z"])
        ligs = lig_map.get(pdb_id, [])
        if not ligs:
            min_dists.append(None)
            druggable.append(0)
            continue

        d = min(_dist((cx, cy, cz), lig) for lig in ligs)
        min_dists.append(float(d))
        druggable.append(1 if d <= p.max_ligand_distance_angstrom else 0)

    df["min_ligand_distance"] = min_dists
    df["druggable"] = druggable
    return df


@dataclass(frozen=True)
class PocketFeaturizationParams:
    neighborhood_radius_angstrom: float = 6.0


def featurize_pockets(
    labeled_pockets: pd.DataFrame,
    protein_ligand_ds: Iterable[dict[str, Any]],
    params: dict[str, Any],
) -> pd.DataFrame:
    """Compute simple chemical-environment features around each pocket centroid.

    Features are computed from protein atoms within a radius of the pocket center.

    Baseline features (per pocket):
      * n_atoms
      * mean_bfactor
      * element counts: C, N, O, S, P, H, other
      * mean distance to centroid

    Added richer features (still fast, heuristic):
      * residue class counts: hydrophobic/polar/pos/neg/aromatic
      * residue counts: n_residues, n_pos_residues, n_neg_residues
      * backbone vs sidechain atom counts
      * distance-weighted element fractions (C/N/O/S)
      * local density proxies (atoms per Å^3, residues per Å^3)
    """
    p = PocketFeaturizationParams(**params)

    # Map pdb_id -> protein structure (Bio.PDB Structure)
    prot_map: dict[str, Any] = {item["pdb_id"]: item["protein"] for item in protein_ligand_ds}

    df = labeled_pockets.copy()

    feats = {
        # existing
        "n_atoms": [],
        "mean_bfactor": [],
        "mean_dist": [],
        "el_C": [],
        "el_N": [],
        "el_O": [],
        "el_S": [],
        "el_P": [],
        "el_H": [],
        "el_other": [],
        # richer
        "n_residues": [],
        "n_backbone_atoms": [],
        "n_sidechain_atoms": [],
        "n_hydrophobic_res": [],
        "n_polar_res": [],
        "n_pos_res": [],
        "n_neg_res": [],
        "n_aromatic_res": [],
        "frac_C_w": [],
        "frac_N_w": [],
        "frac_O_w": [],
        "frac_S_w": [],
        "atom_density": [],
        "residue_density": [],
    }

    # Residue class sets (rough, but useful)
    hydrophobic = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO", "GLY"}
    polar = {"SER", "THR", "ASN", "GLN", "CYS"}
    pos = {"LYS", "ARG", "HIS"}
    neg = {"ASP", "GLU"}
    aromatic = {"PHE", "TYR", "TRP"}
    backbone_names = {"N", "CA", "C", "O", "OXT"}

    # Pre-compute neighborhood volume (sphere) for density
    r = float(p.neighborhood_radius_angstrom)
    sphere_vol = (4.0 / 3.0) * float(np.pi) * (r**3)

    for _, row in df.iterrows():
        pdb_id = row["pdb_id"]
        s = prot_map.get(pdb_id)
        if s is None:
            raise KeyError(f"Protein structure for pdb_id '{pdb_id}' not found")

        center = (float(row["center_x"]), float(row["center_y"]), float(row["center_z"]))
        atoms = list(s.get_atoms())

        coords = np.array([a.coord for a in atoms], dtype=float)
        dists = np.linalg.norm(coords - np.array(center)[None, :], axis=1)
        mask = dists <= p.neighborhood_radius_angstrom

        near_atoms = [a for a, keep in zip(atoms, mask, strict=False) if keep]
        near_dists = dists[mask]

        feats["n_atoms"].append(int(len(near_atoms)))
        feats["mean_dist"].append(float(np.mean(near_dists)) if len(near_dists) else float("nan"))

        bf = [float(getattr(a, "bfactor", 0.0)) for a in near_atoms]
        feats["mean_bfactor"].append(float(np.mean(bf)) if bf else float("nan"))

        # Element counts + distance-weighted fractions
        el_counts = {"C": 0, "N": 0, "O": 0, "S": 0, "P": 0, "H": 0, "other": 0}
        w_counts = {"C": 0.0, "N": 0.0, "O": 0.0, "S": 0.0, "total": 0.0}

        for a in near_atoms:
            el = (getattr(a, "element", "") or "").strip().upper()
            if not el:
                # fallback: guess from atom name's first char
                name = getattr(a, "name", "").strip().upper()
                el = name[:1] if name else ""

            if el in el_counts:
                el_counts[el] += 1
            elif el:
                el_counts["other"] += 1
            else:
                el_counts["other"] += 1

        # Distance weights (closer atoms count more). Add 1.0 to avoid div-by-zero.
        if len(near_atoms):
            weights = 1.0 / (near_dists + 1.0)
            for a, w in zip(near_atoms, weights, strict=False):
                el = (getattr(a, "element", "") or "").strip().upper()
                if not el:
                    name = getattr(a, "name", "").strip().upper()
                    el = name[:1] if name else ""
                if el in {"C", "N", "O", "S"}:
                    w_counts[el] += float(w)
                w_counts["total"] += float(w)

        feats["el_C"].append(el_counts["C"])
        feats["el_N"].append(el_counts["N"])
        feats["el_O"].append(el_counts["O"])
        feats["el_S"].append(el_counts["S"])
        feats["el_P"].append(el_counts["P"])
        feats["el_H"].append(el_counts["H"])
        feats["el_other"].append(el_counts["other"])

        total_w = w_counts["total"]
        if total_w > 0:
            feats["frac_C_w"].append(w_counts["C"] / total_w)
            feats["frac_N_w"].append(w_counts["N"] / total_w)
            feats["frac_O_w"].append(w_counts["O"] / total_w)
            feats["frac_S_w"].append(w_counts["S"] / total_w)
        else:
            feats["frac_C_w"].append(float("nan"))
            feats["frac_N_w"].append(float("nan"))
            feats["frac_O_w"].append(float("nan"))
            feats["frac_S_w"].append(float("nan"))

        # Residue-level features
        near_res_keys: set[tuple[Any, Any, Any]] = set()
        n_backbone = 0
        n_sidechain = 0
        for a in near_atoms:
            res = a.get_parent()
            # ignore waters & non-protein residues for residue-class stats
            if getattr(res, "id", ("", None, ""))[0].strip():
                continue
            resname = getattr(res, "resname", "").strip().upper()
            if not resname:
                continue
            near_res_keys.add((getattr(res, "full_id", None), resname, getattr(res, "id", None)))

            aname = getattr(a, "name", "").strip().upper()
            if aname in backbone_names:
                n_backbone += 1
            else:
                n_sidechain += 1

        # Deduplicate residue names for class counts by residue instance key
        resnames = [rk[1] for rk in near_res_keys]
        feats["n_residues"].append(int(len(resnames)))
        feats["n_backbone_atoms"].append(int(n_backbone))
        feats["n_sidechain_atoms"].append(int(n_sidechain))

        feats["n_hydrophobic_res"].append(int(sum(rn in hydrophobic for rn in resnames)))
        feats["n_polar_res"].append(int(sum(rn in polar for rn in resnames)))
        feats["n_pos_res"].append(int(sum(rn in pos for rn in resnames)))
        feats["n_neg_res"].append(int(sum(rn in neg for rn in resnames)))
        feats["n_aromatic_res"].append(int(sum(rn in aromatic for rn in resnames)))

        # Density proxies
        feats["atom_density"].append(float(len(near_atoms)) / sphere_vol if sphere_vol > 0 else float("nan"))
        feats["residue_density"].append(float(len(resnames)) / sphere_vol if sphere_vol > 0 else float("nan"))

    for k, v in feats.items():
        df[k] = v

    # Make sure label column exists
    if "druggable" not in df.columns:
        raise ValueError("Expected 'druggable' label column in labeled_pockets")

    return df


@dataclass(frozen=True)
class PocketModelParams:
    test_size: float = 0.2
    random_state: int = 42
    xgb: dict[str, Any] | None = None


def split_train_test(
    pocket_features: pd.DataFrame,
    params: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a pocket feature table into train/test tables.

    Keeps all original columns, including `druggable`.

    Params expected (from `params:pocket_model`):
      * test_size
      * random_state
    """
    from sklearn.model_selection import train_test_split

    p = PocketModelParams(**params)

    df = pocket_features.copy()
    if "druggable" not in df.columns:
        raise ValueError("Expected 'druggable' column in pocket_features")

    y = df["druggable"].to_numpy(dtype=int)

    train_df, test_df = train_test_split(
        df,
        test_size=p.test_size,
        random_state=p.random_state,
        stratify=y if len(np.unique(y)) > 1 else None,
    )

    return train_df.reset_index(drop=True), test_df.reset_index(drop=True)


def train_xgb_classifier(
    train_table: pd.DataFrame,
    params: dict[str, Any],
) -> Any:
    """Train a baseline XGBoost classifier to predict `druggable` using the train split."""
    try:
        from xgboost import XGBClassifier
    except Exception as e:  # pragma: no cover
        raise ImportError("xgboost is required for training. Install it with `uv add xgboost`.") from e

    p = PocketModelParams(**params)

    feature_cols = _select_feature_cols(train_table)
    X_train = (
        train_table[feature_cols]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    y_train = train_table["druggable"].to_numpy(dtype=int)

    xgb_params = {
        "n_estimators": 300,
        "max_depth": 5,
        "learning_rate": 0.05,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "reg_lambda": 1.0,
        "objective": "binary:logistic",
        "eval_metric": "logloss",
        "n_jobs": 4,
        "random_state": p.random_state,
    }
    if p.xgb:
        xgb_params.update(p.xgb)

    model = XGBClassifier(**xgb_params)
    model.fit(X_train, y_train)
    return model


def evaluate_classifier(
    model: Any,
    test_table: pd.DataFrame,
    params: dict[str, Any],
) -> tuple[dict[str, Any], Any]:
    """Evaluate model on test split and return metrics + ROC curve figure."""
    from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_fscore_support, roc_curve

    import matplotlib.pyplot as plt

    _ = PocketModelParams(**params)  # currently unused here, but keeps config contract consistent

    feature_cols = _select_feature_cols(test_table)
    X_test = (
        test_table[feature_cols]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    y_test = test_table["druggable"].to_numpy(dtype=int)

    if hasattr(model, "predict_proba"):
        proba = model.predict_proba(X_test)[:, 1]
    else:  # pragma: no cover
        # very defensive fallback
        proba = model.predict(X_test)
        proba = np.asarray(proba, dtype=float)

    pred = (proba >= 0.5).astype(int)

    acc = float(accuracy_score(y_test, pred))
    try:
        auc = float(roc_auc_score(y_test, proba))
    except ValueError:
        auc = float("nan")

    pr, rc, f1, _ = precision_recall_fscore_support(y_test, pred, average="binary", zero_division=0)

    metrics: dict[str, Any] = {
        "n_test": int(len(test_table)),
        "n_test_positive": int(np.sum(y_test)),
        "n_test_negative": int(len(y_test) - np.sum(y_test)),
        "feature_cols": feature_cols,
        "accuracy": acc,
        "roc_auc": auc,
        "precision": float(pr),
        "recall": float(rc),
        "f1": float(f1),
    }

    # Validate JSON-serializability
    json.loads(json.dumps(metrics))

    # ROC curve
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    if np.unique(y_test).size >= 2:
        fpr, tpr, _thr = roc_curve(y_test, proba)
        ax.plot(fpr, tpr, label=f"ROC AUC = {auc:.3f}")
    else:
        ax.text(0.5, 0.5, "ROC undefined (single class in test)", ha="center", va="center")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Pocket druggability ROC curve")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    return metrics, fig


def _select_feature_cols(df: pd.DataFrame) -> list[str]:
    if "druggable" not in df.columns:
        raise ValueError("Expected 'druggable' column")

    drop_cols = {
        "pdb_id",
        "pocket_id",
        "raw",
        "ligand_centroids",
        "min_ligand_distance",
    }
    feature_cols = [c for c in df.columns if c not in drop_cols and c != "druggable"]
    # Keep only numerics
    numeric = [c for c in feature_cols if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric:
        raise ValueError("No numeric feature columns found to train on")
    return numeric


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    ax, ay, az = a
    bx, by, bz = b
    return float(((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2) ** 0.5)
