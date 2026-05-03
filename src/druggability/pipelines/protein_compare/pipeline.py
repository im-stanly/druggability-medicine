from kedro.pipeline import Node, Pipeline  # noqa
from .nodes import find_matching_chains, align_and_calculate_rmsd, predict_binding_pockets, compare_af_pdb_pockets


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([Node(
        func=find_matching_chains,
        inputs=["protein_structure_dataset", "params:aligner_cfg", "params:match_threshold", "decompression_complete"],
        outputs="matched_protein_pairs_dataset",
        name="find_matching_chains_node",
    ),
    Node(
        func=align_and_calculate_rmsd,
        inputs="matched_protein_pairs_dataset",
        outputs=["aligned_protein_pairs_dataset", "protein_alignment_results_dataset", "aligned_completed"],
        name="align_and_calculate_rmsd_node",
    ),
    Node(
        func=predict_binding_pockets,
        inputs=["aligned_paths", "params:p2rank", "params:pocket_probability_threshold", "aligned_completed"],
        outputs="predicted_pockets_dataset",
        name="p2rank_pocket_finding_node",
    ),
    Node(
        func=compare_af_pdb_pockets,
        inputs=["predicted_pockets_dataset", "aligned_protein_pairs_dataset", "params:pocket_match_threshold"], 
        outputs=["pocket_comparison_metrics_dataset", "pocket_comparison_summary_dataset"],
        name="compare_pockets_node",
    )
])
