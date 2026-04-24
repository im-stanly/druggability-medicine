"""
This is a boilerplate pipeline 'protein_parsing'
generated using Kedro 1.3.0
"""

from kedro.pipeline import node, Pipeline  # noqa
from .nodes import do_whatever

def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
        node(
            func=do_whatever,
            inputs="protein_structure_dataset",
            outputs=None
        )
    ])
