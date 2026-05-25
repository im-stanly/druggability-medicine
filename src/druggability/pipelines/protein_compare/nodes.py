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
    for protein_pair in protein_pairs:
        name = protein_pair["name"]
        chain_data_1 = protein_pair["pdb"]
        chain_data_2 = protein_pair["af"]

        chain_1_atoms = [atom for atom in chain_data_1.get_atoms() if atom.get_name() == 'CA']
        chain_2_atoms = [atom for atom in chain_data_2.get_atoms() if atom.get_name() == 'CA']

        superimposer = Superimposer()
        superimposer.set_atoms(chain_1_atoms, chain_2_atoms)
        superimposer.apply(chain_data_1.get_atoms())
        rmses[name] = superimposer.rms
    result_df = pd.DataFrame.from_dict(rmses, orient='index', columns=["rmsd"])
    result_df.reset_index(inplace=True, names=["protein"])

    return result_df


def find_matching_chains(protein_pairs, aligner_cfg, thresh):
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
