"""
Shared constants for the pocket_mining module.

"""

# standard amino acids 

STANDARD_AAS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
})

# residues we skip when parsing (solvent, buffer, crystallisation )
IGNORE_RESIDUES = frozenset({
    "HOH", "DOD", "WAT", "NAG", "MAN", "UNK", "GLC", "ABA",
    "MPD", "GOL", "SO4", "PO4", "EDO", "ACT", "PEG", "BME",
    "DMS", "FMT", "EPE", "CIT", "TRS",
})

# residue groups (by side-chain chemistry)

HYDROPHOBIC = frozenset({"ALA","VAL","LEU","ILE","PRO","PHE","TRP","MET","CYS"})
POLAR       = frozenset({"SER","THR","ASN","GLN","TYR","HIS"})
CHARGED_POS = frozenset({"LYS","ARG"})
CHARGED_NEG = frozenset({"ASP","GLU"})
AROMATIC    = frozenset({"PHE","TYR","TRP","HIS"})

# physicochemical scales

HYDROPHOBICITY = {
    "ALA": 0.62, "ARG":-2.53, "ASN":-0.78, "ASP":-0.90,
    "CYS": 0.29, "GLN":-0.85, "GLU":-0.74, "GLY": 0.48,
    "HIS":-0.40, "ILE": 1.38, "LEU": 1.06, "LYS":-1.50,
    "MET": 0.64, "PHE": 1.19, "PRO": 0.12, "SER":-0.18,
    "THR":-0.05, "TRP": 0.81, "TYR": 0.26, "VAL": 1.08,
}

RESIDUE_VOLUME = {
    "ALA": 88.6, "ARG":173.4, "ASN":114.1, "ASP":111.1,
    "CYS":108.5, "GLN":143.8, "GLU":138.4, "GLY": 60.1,
    "HIS":153.2, "ILE":166.7, "LEU":166.7, "LYS":168.6,
    "MET":162.9, "PHE":189.9, "PRO":112.7, "SER": 89.0,
    "THR":116.1, "TRP":227.8, "TYR":193.6, "VAL":140.0,
}

VDW_RADII = {
    "H":1.20,"C":1.70,"N":1.55,"O":1.52,
    "S":1.80,"P":1.80,"F":1.47,"CL":1.75,
    "MG":1.73,"CA":2.31,"ZN":1.39,"FE":1.56,
    "MN":1.61,"NA":2.27,"K":2.75,"SE":1.90,
}

# params

PROBE_RADIUS   = 1.4   # water probe for solvent-accessible surface (Å)
DEFAULT_N_POINTS = 2000  # surface points per protein
FEATURE_RADIUS = 8.0   # neighbourhood for feature lookups (Å)

POCKET_RADIUS     = 5.0   # point is "pocket" if ≤ this from any ligand atom
NON_POCKET_RADIUS = 10.0  # point is "not pocket" if ≥ this from all ligands

PROBABILITY_THRESHOLD = 0.5
CLUSTER_RADIUS        = 5.0
MIN_CLUSTER_SIZE      = 10
