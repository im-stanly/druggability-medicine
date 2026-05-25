"""Features that describe a spot on the protein surface.

For each surface point we look at all protein atoms within 8 Å and
summarise what we see into 28 numbers. A Random Forest uses these
to decide: pocket or not pocket?
"""

import numpy as np
from scipy.spatial import KDTree

from .constants import (
    FEATURE_RADIUS,
    HYDROPHOBIC, POLAR, CHARGED_POS, CHARGED_NEG, AROMATIC,
    HYDROPHOBICITY, RESIDUE_VOLUME,
)
from .parser import ProteinAtoms


def featurize(
    surface_points: np.ndarray,
    protein: ProteinAtoms,
    radius: float = FEATURE_RADIUS,
) -> np.ndarray:
    """Turn surface points into a (N, 28) feature matrix.

    The 28 features:
      [ 0– 4]  atom type fractions: C, N, O, S, other
      [ 5–10]  residue class fractions: hydrophobic, polar, pos, neg, aromatic, gly
      [11]      avg hydrophobicity (Eisenberg)
      [12]      avg residue volume (Å³)
      [13]      atom density (atoms / Å³)
      [14]      fraction of surface-exposed atoms
      [15]      distance to nearest CA (Å)
      [16]      protrusion (negative = dent, positive = bump)
      [17]      fraction backbone atoms (N, CA, C, O)
      [18]      fraction sidechain atoms
      [19–21]   mean offset of nearby atoms (x, y, z)
      [22–24]   spatial spread of nearby atoms (x, y, z)  ← top features
      [25]      log(number of nearby atoms)
      [26]      fraction of atoms from charged residues
      [27]      net charge balance (pos − neg)
    """
    tree = KDTree(protein.coords)
    n_points = len(surface_points)
    backbone = {"N", "CA", "C", "O"}
    sphere_vol = (4.0 / 3.0) * np.pi * radius ** 3

    # --- precompute per-atom masks (once) ---
    elem = np.array(protein.elements)
    resn = np.array(protein.residue_names)

    is_C = elem == "C"
    is_N = elem == "N"
    is_O = elem == "O"
    is_S = elem == "S"

    is_hyd  = np.array([r in HYDROPHOBIC for r in resn])
    is_pol  = np.array([r in POLAR for r in resn])
    is_pos  = np.array([r in CHARGED_POS for r in resn])
    is_neg  = np.array([r in CHARGED_NEG for r in resn])
    is_aro  = np.array([r in AROMATIC for r in resn])
    is_gly  = resn == "GLY"

    hydro  = np.array([HYDROPHOBICITY.get(r, 0.0) for r in resn])
    volume = np.array([RESIDUE_VOLUME.get(r, 0.0) for r in resn])

    # protein centre + mean radius (for protrusion / surface exposure)
    centre = protein.coords.mean(axis=0)
    mean_r = np.linalg.norm(protein.coords - centre, axis=1).mean()

    # CA atoms for distance feature
    ca_mask = np.array([n == "CA" for n in protein.atom_names])
    ca_tree = KDTree(protein.coords[ca_mask]) if ca_mask.any() else None

    features = np.zeros((n_points, 28))

    for i, pt in enumerate(surface_points):
        idx = tree.query_ball_point(pt, radius)
        n = len(idx)
        if n == 0:
            features[i, 15] = radius   # far from CA
            continue

        local = protein.coords[idx]
        rel = local - pt
        names = [protein.atom_names[j] for j in idx]
        dists = np.linalg.norm(local - centre, axis=1)

        # 0–4: atom type fractions
        n_C = is_C[idx].sum()
        n_N = is_N[idx].sum()
        n_O = is_O[idx].sum()
        n_S = is_S[idx].sum()
        features[i, 0] = n_C / n
        features[i, 1] = n_N / n
        features[i, 2] = n_O / n
        features[i, 3] = n_S / n
        features[i, 4] = (n - n_C - n_N - n_O - n_S) / n

        # 5–10: residue class fractions
        features[i, 5]  = is_hyd[idx].sum() / n
        features[i, 6]  = is_pol[idx].sum() / n
        features[i, 7]  = is_pos[idx].sum() / n
        features[i, 8]  = is_neg[idx].sum() / n
        features[i, 9]  = is_aro[idx].sum() / n
        features[i, 10] = is_gly[idx].sum() / n

        # 11–12: physicochemical averages
        features[i, 11] = hydro[idx].mean()
        features[i, 12] = volume[idx].mean()

        # 13: atom density
        features[i, 13] = n / sphere_vol

        # 14: surface-exposed fraction
        features[i, 14] = (dists > mean_r).sum() / n

        # 15: distance to nearest CA
        features[i, 15] = ca_tree.query(pt, k=1)[0] if ca_tree else radius

        # 16: protrusion
        features[i, 16] = dists.mean() - mean_r

        # 17–18: backbone / sidechain
        n_bb = sum(1 for nm in names if nm.strip() in backbone)
        features[i, 17] = n_bb / n
        features[i, 18] = (n - n_bb) / n

        # 19–21: mean offset of neighbourhood
        features[i, 19:22] = rel.mean(axis=0)

        # 22–24: spatial spread (top features by importance)
        if n > 1:
            features[i, 22:25] = rel.std(axis=0)

        # 25: log atom count
        features[i, 25] = np.log(n)

        # 26–27: charge
        pos_n = is_pos[idx].sum()
        neg_n = is_neg[idx].sum()
        features[i, 26] = (pos_n + neg_n) / n
        features[i, 27] = (pos_n - neg_n) / n

    return features


def feature_names() -> list[str]:
    return [
        "C_frac", "N_frac", "O_frac", "S_frac", "other_frac",
        "hydrophobic", "polar", "positive", "negative",
        "aromatic", "glycine",
        "avg_hydrophobicity", "avg_volume",
        "atom_density", "surf_exposed",
        "dist_to_CA", "protrusion",
        "frac_backbone", "frac_sidechain",
        "offset_x", "offset_y", "offset_z",
        "spread_x", "spread_y", "spread_z",
        "log_n_atoms", "frac_charged", "net_charge",
    ]
