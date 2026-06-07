import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import math

def compare_af_pdb_pockets(pockets_df: pd.DataFrame, protein_dataset, pocket_match_threshold: float, pocket_shared_residues_number: int) -> pd.DataFrame:
    if pockets_df.empty:
        return pd.DataFrame(), pd.DataFrame() # Poprawka: Kedro oczekuje dwóch outputów wg sygnatury w return

    comparison_results = []

    for protein_name, group in pockets_df.groupby("protein_name"):
        pdb_pockets = group[group["source_type"] == "PDB"].copy()
        af_pockets = group[group["source_type"] == "AF"].copy()
        
        if pdb_pockets.empty or af_pockets.empty:
            continue

        pdb_centers = pdb_pockets[['center_x', 'center_y', 'center_z']].values
        af_centers = af_pockets[['center_x', 'center_y', 'center_z']].values
        dist_matrix = cdist(pdb_centers, af_centers, metric='euclidean')
        
        pair_data = next((p for p in protein_dataset if p["name"] == protein_name), None)

        for i, pdb_row in enumerate(pdb_pockets.itertuples()):
            closest_af_idx = np.argmin(dist_matrix[i])
            center_distance = dist_matrix[i, closest_af_idx]
            
            if center_distance <= pocket_match_threshold:
                af_row = af_pockets.iloc[closest_af_idx]
                
                pocket_rmsd = None
                avg_plddt = None
                
                if pair_data:
                    pdb_residue_strings = str(pdb_row.residue_ids).split()
                    af_residue_strings = str(af_row.residue_ids).split()
                    shared_residues = set(pdb_residue_strings).intersection(set(af_residue_strings))
                    
                    sq_distances = []
                    plddt_scores = [] 
                    
                    if len(shared_residues) < pocket_shared_residues_number:
                        comparison_results.append({
                            "protein_name": protein_name,
                            "pdb_pocket_rank": pdb_row.rank,
                            "af_pocket_rank": af_row["rank"],
                            "pdb_probability": pdb_row.probability,
                            "af_probability": af_row["probability"],
                            "center_shift_A": center_distance,
                            "pocket_backbone_rmsd_A": None,
                            "avg_pocket_plddt": None, 
                            "shared_residue_count": len(shared_residues),
                            "match_status": f"Insufficient shared residues ({len(shared_residues)})"
                        })
                        continue

                    for res_str in shared_residues:
                        try:
                            chain_id, res_id = res_str.split("_")
                            res_id = int(res_id)
                            pdb_ca = pair_data["pdb"][0][chain_id][res_id]['CA']
                            af_ca = pair_data["af"][0][chain_id][res_id]['CA']
                            
                            diff = pdb_ca.coord - af_ca.coord
                            sq_distances.append(np.sum(diff ** 2))
                            
                            # getting pLDDT (saved in AlphaFold's model B-factor)
                            plddt_scores.append(af_ca.get_bfactor())
                            
                        except (KeyError, ValueError):
                            continue 
                    
                    if sq_distances:
                        pocket_rmsd = math.sqrt(sum(sq_distances) / len(sq_distances))
                    if plddt_scores:
                        avg_plddt = sum(plddt_scores) / len(plddt_scores)

                comparison_results.append({
                    "protein_name": protein_name,
                    "pdb_pocket_rank": pdb_row.rank,
                    "af_pocket_rank": af_row["rank"],
                    "pdb_probability": pdb_row.probability,
                    "af_probability": af_row["probability"],
                    "center_shift_A": center_distance,
                    "pocket_backbone_rmsd_A": pocket_rmsd,
                    "avg_pocket_plddt": avg_plddt,
                    "shared_residue_count": len(shared_residues) if pair_data else 0,
                    "match_status": "Matched"
                })

            else:
                comparison_results.append({
                    "protein_name": protein_name,
                    "pdb_pocket_rank": pdb_row.rank,
                    "af_pocket_rank": None,
                    "pdb_probability": pdb_row.probability,
                    "af_probability": None,
                    "center_shift_A": None,
                    "pocket_backbone_rmsd_A": None,
                    "avg_pocket_plddt": None, 
                    "shared_residue_count": 0,
                    "match_status": "Missing in AF"
                })

    comparison_df = pd.DataFrame(comparison_results)

    if comparison_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    matched_pockets = comparison_df[comparison_df["match_status"] == "Matched"]
    avg_df = matched_pockets.groupby("protein_name").agg(
        avg_pocket_rmsd_A=("pocket_backbone_rmsd_A", "mean"),
        avg_center_shift_A=("center_shift_A", "mean"),
        avg_protein_pocket_plddt=("avg_pocket_plddt", "mean"), # NOWE
        total_pockets=("protein_name", "size")
    ).reset_index()

    return comparison_df, avg_df

from scipy.stats import spearmanr, pearsonr

def verify_hypothesis(comparison_df: pd.DataFrame) -> dict:
    """
    Verify the following hypotheses based on the pocket comparison results:
    1. >80% matches have RMSD < 2 Å.
    2. RMSD is lower for the pockets with higher pLDDT.
    """
    matched = comparison_df[
        (comparison_df["match_status"] == "Matched") & 
        (comparison_df["pocket_backbone_rmsd_A"].notnull()) &
        (comparison_df["avg_pocket_plddt"].notnull())
    ]
    
    if matched.empty:
        return {"error": "No matching pockets found."}

    total_matched = len(matched)
    rmsd_under_2 = len(matched[matched["pocket_backbone_rmsd_A"] < 2.0])
    percent_under_2 = (rmsd_under_2 / total_matched) * 100

    # --- PART 2: Decline in accuracy for low pLDDT ---
    # Correlation between pLDDT and RMSD. We expect a negative correlation 
    # (the higher the pLDDT, the lower the RMSD/error).
    plddt_values = matched["avg_pocket_plddt"]
    rmsd_values = matched["pocket_backbone_rmsd_A"]
    
    pearson_corr, p_value = pearsonr(plddt_values, rmsd_values)
    spearman_corr, s_p_value = spearmanr(plddt_values, rmsd_values)
    
    # Category for better visualization (binning pLDDT)
    matched_copy = matched.copy()
    matched_copy["plddt_bin"] = pd.cut(
        matched_copy["avg_pocket_plddt"], 
        bins=[0, 50, 70, 90, 100], 
        labels=["Very low (<50)", "Low (50-70)", "Good (70-90)", "Very good (>90)"]
    )
    rmsd_by_bin = matched_copy.groupby("plddt_bin")["pocket_backbone_rmsd_A"].mean().to_dict()

    return {
        "Hypothesis 1 (>80% have RMSD < 2A)": {
            "Total matched": total_matched,
            "RMSD below 2": rmsd_under_2,
            "Percentage": round(percent_under_2, 2)
        },
        "Hypothesis 2 (Correlation between pLDDT and RMSD)": {
            "Pearson correlation": round(pearson_corr, 3),
            "Pearson p-value": p_value,
            "Spearman correlation": round(spearman_corr, 3),
            "Average RMSD by pLDDT": rmsd_by_bin
        }
    }