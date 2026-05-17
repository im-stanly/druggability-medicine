"""Pipeline that builds a baseline pocket druggability model.

Stages:
  1) label pockets by distance to nearest ligand centroid
  2) featurize pocket centroid neighborhood from protein structure
  3) split into train/test
  4) train an XGBoost classifier
  5) evaluate + save ROC curve
"""

from kedro.pipeline import Node, Pipeline

from .nodes import (
    evaluate_classifier,
    featurize_pockets,
    label_pockets_by_ligand_distance,
    split_train_test,
    train_xgb_classifier,
)


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline(
        [
            Node(
                func=label_pockets_by_ligand_distance,
                inputs=dict(
                    pockets="p2rank_pockets",
                    protein_ligand_ds="protein_ligand_dataset",
                    params="params:pocket_labeling",
                ),
                outputs="labeled_pockets",
                name="label_pockets_by_ligand_distance",
            ),
            Node(
                func=featurize_pockets,
                inputs=dict(
                    labeled_pockets="labeled_pockets",
                    protein_ligand_ds="protein_ligand_dataset",
                    params="params:pocket_featurization",
                ),
                outputs="pocket_features",
                name="featurize_pockets",
            ),
            Node(
                func=split_train_test,
                inputs=dict(
                    pocket_features="pocket_features",
                    params="params:pocket_model",
                ),
                outputs=["pocket_train_table", "pocket_test_table"],
                name="split_train_test",
            ),
            Node(
                func=train_xgb_classifier,
                inputs=dict(
                    train_table="pocket_train_table",
                    params="params:pocket_model",
                ),
                outputs="pocket_xgb_model",
                name="train_xgb_classifier",
            ),
            Node(
                func=evaluate_classifier,
                inputs=dict(
                    model="pocket_xgb_model",
                    test_table="pocket_test_table",
                    params="params:pocket_model",
                ),
                outputs=["pocket_model_metrics", "pocket_roc_curve"],
                name="evaluate_classifier",
            ),
        ]
    )
