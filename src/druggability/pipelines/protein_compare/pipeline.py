from kedro.pipeline import Node, Pipeline  # noqa
from .nodes import find_matching_chains, align_and_calculate_rmsd


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([Node(
        func=find_matching_chains,
        inputs=["protein_structure_dataset", "params:aligner_cfg", "params:match_threshold"],
        outputs="matched_protein_pairs_dataset",
        name="find_matching_chains_node",
    ),
    Node(
        func=align_and_calculate_rmsd,
        inputs="matched_protein_pairs_dataset",
        outputs="protein_alignment_results_dataset",
        name="align_and_calculate_rmsd_node",
    )
    ])
