"""Protein surface point generation.

Scatters ~2000 points on the solvent-accessible surface, then labels
each one as pocket / non-pocket / ambiguous based on ligand distance.
"""

import numpy as np
from scipy.spatial import KDTree

from .constants import PROBE_RADIUS, DEFAULT_N_POINTS, VDW_RADII, POCKET_RADIUS, NON_POCKET_RADIUS
from .parser import ProteinAtoms


def _surface_atoms(protein: ProteinAtoms, probe: float = PROBE_RADIUS) -> np.ndarray:
    """Return indices of surface-exposed atoms.

    An atom is "surface" if the number of neighbouring atoms within
    vdW+probe distance is below a threshold.  This is a fast geometric
    proxy for SASA — more neighbours = more buried.
    """
    coords = protein.coords
    n = len(coords)
    tree = KDTree(coords)

    # Count neighbours within vdW + probe for each atom
    neighbour_counts = np.zeros(n, dtype=int)
    for i in range(n):
        r_i = VDW_RADII.get(protein.elements[i].upper(), 1.70) + probe
        neighbour_counts[i] = len(tree.query_ball_point(coords[i], r_i)) - 1

    # An atom is a "surface" if it has ≤ 6 neighbours at vdW+probe distance.
    # (Buried atoms typically have 8-12 neighbours; surface atoms have 2-6.)
    # Also always include CA atoms since they trace the backbone surface.
    is_ca = np.array([name == "CA" for name in protein.atom_names])
    surface = (neighbour_counts <= 6) | is_ca
    return np.where(surface)[0]


def generate_surface_points(
    protein: ProteinAtoms,
    n_points: int = DEFAULT_N_POINTS,
    probe_radius: float = PROBE_RADIUS,
    random_seed: int | None = 42,
) -> np.ndarray:
    """Generate points on the solvent-accessible surface.

    1. Find surface-exposed atoms via local neighbour count.
    2. Scatter candidate points at vdW+probe distance in random directions.
    3. Reject points that clash with any protein atom.
    4. Subsample to n_points.
    """
    rng = np.random.RandomState(random_seed)
    coords = protein.coords

    surface_idx = _surface_atoms(protein, probe_radius)
    if len(surface_idx) == 0:
        surface_idx = np.arange(len(coords))
    tree = KDTree(coords)
    n_per_atom = max(2, 10 * n_points // len(surface_idx))
    candidates = []

    for idx in surface_idx:
        r_vdw = VDW_RADII.get(protein.elements[idx].upper(), 1.70)
        center = coords[idx]

        for _ in range(n_per_atom):
            d = rng.randn(3)
            d /= np.linalg.norm(d)
            pt = center + d * (r_vdw + probe_radius)

            dist, nearest = tree.query(pt, k=1)
            r_near = VDW_RADII.get(protein.elements[nearest].upper(), 1.70)
            if dist >= r_near + probe_radius:
                candidates.append(pt)

    if not candidates:
        candidates = [coords[i] + rng.randn(3) * 2.5 for i in surface_idx]

    candidates = np.array(candidates)
    if len(candidates) > n_points:
        candidates = candidates[rng.choice(len(candidates), n_points, replace=False)]
    elif len(candidates) < n_points and len(candidates) > 0:
        candidates = candidates[rng.choice(len(candidates), n_points, replace=True)]

    return candidates


def label_points(
    surface_points: np.ndarray,
    ligand_coords: np.ndarray,
    pocket_radius: float = POCKET_RADIUS,
    non_pocket_radius: float = NON_POCKET_RADIUS,
) -> np.ndarray:
    """Label points: 1 = pocket, 0 = non-pocket, -1 = ambiguous (ignore)."""
    if len(ligand_coords) == 0:
        return np.full(len(surface_points), -1, dtype=int)

    dists, _ = KDTree(ligand_coords).query(surface_points, k=1)
    labels = np.full(len(surface_points), -1, dtype=int)
    labels[dists <= pocket_radius] = 1
    labels[dists >= non_pocket_radius] = 0
    return labels
