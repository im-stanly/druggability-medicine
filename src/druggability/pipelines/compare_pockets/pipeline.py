from kedro.pipeline import Node, Pipeline  # noqa

from .nodes import compare_af_pdb_pockets, verify_hypothesis


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
    Node(
        func=compare_af_pdb_pockets,
        inputs=["predicted_pockets_dataset", "aligned_protein_pairs_dataset", "params:pocket_match_threshold"], 
        outputs=["pocket_comparison_metrics_dataset", "pocket_comparison_summary_dataset"],
        name="compare_pockets_node",
    ),
    Node(
        func=verify_hypothesis,
        inputs="pocket_comparison_metrics_dataset",
        outputs="hypothesis_verification_results",
        name="hypothesis_verification_node",
    )])
