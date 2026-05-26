"""Project pipelines."""

from __future__ import annotations

from druggability.pipelines.protein_unzip import create_pipeline as create_protein_unzip_pipeline
from druggability.pipelines.protein_compare import create_pipeline as create_protein_compare_pipeline
from kedro.pipeline import Pipeline
from druggability.pipelines.extract_pockets.pipeline import create_pipeline as create_extract_pockets_pipeline
from druggability.pipelines.pocket_ml.pipeline import create_pipeline as create_pocket_ml_pipeline


def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    unzip_pipeline = create_protein_unzip_pipeline()
    compare_pipeline = create_protein_compare_pipeline()
    extract_pockets_pipeline = create_extract_pockets_pipeline()
    model_pipeline = create_pocket_ml_pipeline()
    train_model_pipeline = extract_pockets_pipeline + model_pipeline

    pipelines = {"__default__": train_model_pipeline,
                 "model_training": train_model_pipeline,
                 "unzip": unzip_pipeline,
                 "compare": compare_pipeline,
                 "extract_pockets": extract_pockets_pipeline,
                 "pocket_ml": model_pipeline
                 }
    return pipelines
