from pathlib import Path
import shutil
import shutil
import subprocess
from sys import platform
from typing import Any, Callable, Dict
import platform as sys_platform

from kedro_datasets import partitions
import pandas as pd
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
from Bio.PDB import Select, Superimposer

def _get_chains(structure):
    chains = {}
    model = structure[0]

    for chain in model:
        residues = []
        sequence = ""
        for residue in chain:
            if is_aa(residue, standard=True):
                try:
                    seq_char = seq1(residue.resname)
                    sequence += seq_char
                    residues.append(residue)
                except KeyError: # Non-standard amino acid, skip
                    continue
        if sequence:
            chains[chain.id] = (residues, sequence)

    return chains

def _get_aligner(mode, open_gap_score, extend_gap_score):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    return aligner

def _find_matching_chains(chain_data_1, chain_data_2, aligner, match_threshold):
    keep_res1, keep_res2 = [], []
    ca_atoms_1, ca_atoms_2 = [], []
    matched_chains_1 = set()
    matched_chains_2 = set()
    for chain_id_1, (residues_1, seq_1) in chain_data_1.items():
        if chain_id_1 in matched_chains_1:
            continue
        for chain_id_2, (residues_2, seq_2) in chain_data_2.items():
            if chain_id_2 in matched_chains_2:
                continue
            alignments = aligner.align(seq_1, seq_2)
            best_alignment = alignments[0]

            max_possible_score = min(len(seq_1), len(seq_2))
            percentage_matched = best_alignment.score / max_possible_score if max_possible_score > 0 else 0

            if percentage_matched >= match_threshold:
                aligned_seq_1, aligned_seq_2 = best_alignment.aligned
                matched_chains_1.add(chain_id_1)
                matched_chains_2.add(chain_id_2)
                for (start1, end1), (start2, end2) in zip(aligned_seq_1, aligned_seq_2):
                    keep_res1.extend(residues_1[start1:end1])
                    keep_res2.extend(residues_2[start2:end2])

    return keep_res1, keep_res2


def align_and_calculate_rmsd(protein_pairs):
    rmses = {}
    aligned_pairs = []
    for protein_pair in protein_pairs:
        name = protein_pair["name"]
        chain_data_1 = protein_pair["pdb"]
        chain_data_2 = protein_pair["af"]

        chain_1_atoms = [atom for atom in chain_data_1.get_atoms() if atom.get_name() == 'CA']
        chain_2_atoms = [atom for atom in chain_data_2.get_atoms() if atom.get_name() == 'CA']

        superimposer = Superimposer()
        superimposer.set_atoms(chain_1_atoms, chain_2_atoms)
        superimposer.apply(chain_data_2.get_atoms())

        aligned_pairs.append({
            "name": name,
            "pdb": chain_data_1,
            'af': chain_data_2
        })

        rmses[name] = superimposer.rms
    result_df = pd.DataFrame.from_dict(rmses, orient='index', columns=["rmsd"])
    result_df.reset_index(inplace=True, names=["protein"])

    return aligned_pairs, result_df, []


def find_matching_chains(protein_pairs, aligner_cfg, thresh, _):
    aligner = _get_aligner(**aligner_cfg)
    matched_pairs = []
    for data_dict in protein_pairs:
        pdb_prot, af_prot = data_dict["pdb"], data_dict["af"]
        chain_data_1 = _get_chains(pdb_prot)
        chain_data_2 = _get_chains(af_prot)
        keep_res1, keep_res2 = _find_matching_chains(chain_data_1, chain_data_2, aligner, thresh)
        if keep_res1 and keep_res2:
            matched_pairs.append({
                "name": data_dict["name"],
                "keep_pdb": keep_res1,
                "pdb": pdb_prot,
                "keep_af": keep_res2,
                "af": af_prot,
            })

    return matched_pairs

def predict_binding_pockets(partitions: Dict[str, Callable[[], Any]], p2rank_params: dict, pocket_probability_threshold: float, _) -> pd.DataFrame:
    """
    Kedro node to run P2Rank on aligned protein structures.
    """
    p2rank_exec = Path(p2rank_params["executable_path"]).resolve()
    work_dir = Path(p2rank_params["working_dir"]).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    dataset_file = work_dir / "batch_run.ds"
    

    if sys_platform.system() == "Windows" and p2rank_exec.suffix != ".bat":
        p2rank_exec = p2rank_exec.with_name(f"{p2rank_exec.name}.bat")

    # We create a specific staging folder so we don't permanently rename your raw data
    staging_dir = work_dir / "staging_cifs"
    output_dir = work_dir / "p2rank_outputs"
    staging_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    # 3. Find all .cif files recursively and stage them
    # rglob("*.cif") searches the base directory and ALL subdirectories automatically
    staging_dir.mkdir(parents=True, exist_ok=True)

    with open(dataset_file, "w") as f:
        for name, path_func in partitions.items():
            path = Path(path_func()).resolve()
            protein_name = path.parent.name 
            source_type = path.stem
            unique_name = f"{protein_name}_{source_type}.cif"
            staged_path = staging_dir / unique_name
            
            # Copy the file to our temporary staging area
            shutil.copy2(path, staged_path)
            f.write(f"{staged_path}\n")

    command = [
        str(p2rank_exec),
        "predict",
        "-o", str(output_dir),    
        str(dataset_file),
    ]

    try:
        subprocess.run(command, capture_output=True, text=True, check=True, timeout=300)
    except subprocess.CalledProcessError as e:
        print(f"P2Rank failed for {dataset_file}:\n{e.stderr}")

    csv_files = list(output_dir.rglob("*_predictions.csv"))
    if not csv_files:
        print(f"No prediction CSVs found in {output_dir}")
        return pd.DataFrame()
        
    all_dfs = []
    for csv_path in csv_files:
        # 2. Read the CSV 
        # P2Rank leaves weird spaces after commas, skipinitialspace fixes this
        df = pd.read_csv(csv_path, skipinitialspace=True)
        
        if df.empty:
            continue
            
        df = df[df["probability"].gt(pocket_probability_threshold)]
        # 3. Parse the filename
        # Example: "BRCA1_mutant_af.cif_predictions.csv" -> "BRCA1_mutant_af"
        clean_name = csv_path.name.replace("_predictions.csv", "").replace(".cif", "")
        
        # Split from the right at the LAST underscore
        parts = clean_name.rsplit("_", 1)
        
        # 4. Inject the metadata columns
        df["protein_name"] = parts[0] if len(parts) > 1 else clean_name
        df["source_type"] = parts[1] if len(parts) > 1 else "unknown"
        
        all_dfs.append(df)

    if all_dfs:
        print(f"Successfully merged {len(csv_files)} files.")
        return pd.concat(all_dfs, ignore_index=True)
    else:
        return pd.DataFrame()

import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import math

def compare_af_pdb_pockets(pockets_df: pd.DataFrame, protein_dataset, pocket_match_threshold: float) -> pd.DataFrame:
    """
    Matches PDB and AF pockets based on spatial proximity and calculates 
    both the Center Shift and the structural RMSD of the pocket residues.
    """
    if pockets_df.empty:
        return pd.DataFrame()

    comparison_results = []

    # Group by protein so we only compare a protein to its own AF model
    for protein_name, group in pockets_df.groupby("protein_name"):
        
        pdb_pockets = group[group["source_type"] == "PDB"].copy()
        af_pockets = group[group["source_type"] == "AF"].copy()
        
        if pdb_pockets.empty or af_pockets.empty:
            continue # Skip if one structure didn't produce any pockets

        # 1. MATCH THE POCKETS using Euclidean Distance Matrix
        # Extract the X, Y, Z coordinates into numpy arrays
        pdb_centers = pdb_pockets[['center_x', 'center_y', 'center_z']].values
        af_centers = af_pockets[['center_x', 'center_y', 'center_z']].values
        
        # Calculate distance between every PDB pocket and every AF pocket
        # dist_matrix[i, j] will be the distance between PDB pocket i and AF pocket j
        dist_matrix = cdist(pdb_centers, af_centers, metric='euclidean')
        
        # Find the Biopython structures for atom-level RMSD
        pair_data = next((p for p in protein_dataset if p["name"] == protein_name), None)

        # 2. EVALUATE MATCHES
        for i, pdb_row in enumerate(pdb_pockets.itertuples()):
            # Find the index of the closest AF pocket
            closest_af_idx = np.argmin(dist_matrix[i])
            center_distance = dist_matrix[i, closest_af_idx]
            
            if center_distance <= pocket_match_threshold:
                af_row = af_pockets.iloc[closest_af_idx]
                
                pocket_rmsd = None
                shared_residues = [] # Initialize to avoid errors
                if pair_data:
                    pdb_residue_strings = str(pdb_row.residue_ids).split()
                    af_residue_strings = str(af_row.residue_ids).split()
                    shared_residues = set(pdb_residue_strings).intersection(set(af_residue_strings))
                    
                    sq_distances = []
                    for res_str in shared_residues:
                        try:
                            chain_id, res_id = res_str.split("_")
                            res_id = int(res_id)
                            pdb_ca = pair_data["pdb"][0][chain_id][res_id]['CA']
                            af_ca = pair_data["af"][0][chain_id][res_id]['CA']
                            diff = pdb_ca.coord - af_ca.coord
                            sq_distances.append(np.sum(diff ** 2))
                        except (KeyError, ValueError):
                            continue 
                    
                    if sq_distances:
                        pocket_rmsd = math.sqrt(sum(sq_distances) / len(sq_distances))

                # Append the matched data
                comparison_results.append({
                    "protein_name": protein_name,
                    "pdb_pocket_rank": pdb_row.rank,
                    "af_pocket_rank": af_row["rank"],
                    "pdb_probability": pdb_row.probability,
                    "af_probability": af_row["probability"],
                    "center_shift_A": center_distance,
                    "pocket_backbone_rmsd_A": pocket_rmsd,
                    "shared_residue_count": len(shared_residues) if pair_data else 0,
                    "match_status": "Matched" # Add a helpful tag!
                })

            # Scenario B: NO MATCH FOUND (AlphaFold missed it)
            else:
                comparison_results.append({
                    "protein_name": protein_name,
                    "pdb_pocket_rank": pdb_row.rank,
                    "af_pocket_rank": None,        # Blank because AF missed it
                    "pdb_probability": pdb_row.probability,
                    "af_probability": None,        # Blank
                    "center_shift_A": None,        # Blank
                    "pocket_backbone_rmsd_A": None,# Blank
                    "shared_residue_count": 0,
                    "match_status": "Missing in AF" # Flag it as a missing pocket
                })
    import pandas as pd
    comparison_df = pd.DataFrame(comparison_results)



    matched_pockets = comparison_df[comparison_df["match_status"] == "Matched"]
    avg_df = matched_pockets.groupby("protein_name").agg(
        avg_pocket_rmsd_A=("pocket_backbone_rmsd_A", "mean"),
        avg_center_shift_A=("center_shift_A", "mean"),
        
        # .size counts ALL rows (total pockets)
        total_pockets=("protein_name", "size")).reset_index()

    return comparison_df, avg_df