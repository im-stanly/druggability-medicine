"""Kedro pipeline for training the binding pocket classifier."""

from kedro.pipeline import Pipeline, node

from .nodes import train_pocket_classifier


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
        node(
            func=train_pocket_classifier,
            inputs=[
                "params:pocket_model.data_dir",
                "params:pocket_model.n_points",
            ],
            outputs=["pocket_model_metrics", "pocket_model"],
            name="train_pocket_classifier_node",
        ),
    ])
