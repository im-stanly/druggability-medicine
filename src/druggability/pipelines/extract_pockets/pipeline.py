"""
This is a boilerplate pipeline 'extract_pockets'
generated using Kedro 1.3.0
"""

from kedro.pipeline import Node, Pipeline  # noqa
from .nodes import run_p2rank_and_parse_pockets


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
        Node(
            func=run_p2rank_and_parse_pockets,
            inputs=dict(protein_ligand_ds="protein_ligand_dataset", p2rank="params:p2rank"),
            outputs="p2rank_pockets",
            name="run_p2rank_and_parse_pockets",
        )
    ])
