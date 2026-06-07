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

    dataset_pdb = work_dir / "batch_pdb.ds"
    dataset_af = work_dir / "batch_af.ds"
    pdb_files = []
    af_files = []

    # 1. Rozdzielenie plików na dwie grupy na podstawie 'stem'
    for name, path_func in partitions.items():
        path = Path(path_func()).resolve()
        protein_name = path.parent.name 
        source_type = path.stem
        
        unique_name = f"{protein_name}_{source_type}.cif"
        staged_path = staging_dir / unique_name
        
        shutil.copy2(path, staged_path)
        
        if source_type.lower().startswith("af"):
            af_files.append(str(staged_path))
        else:
            pdb_files.append(str(staged_path))

    def run_p2rank(dataset_path, file_paths, extra_args):
        if not file_paths:
            return 
        
        with open(dataset_path, "w") as f:
            f.write("\n".join(file_paths) + "\n")

        command = [str(p2rank_exec), "predict"] + extra_args + ["-o", str(output_dir), str(dataset_path)]
        
        try:
            subprocess.run(command, capture_output=True, text=True, check=True, timeout=300)
            print(f"P2Rank success for {dataset_path.name}")
        except subprocess.CalledProcessError as e:
            print(f"P2Rank failed for {dataset_path.name}:\n{e.stderr}")

    run_p2rank(dataset_pdb, pdb_files, []) 
    run_p2rank(dataset_af, af_files, ["-c", "alphafold"])

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

