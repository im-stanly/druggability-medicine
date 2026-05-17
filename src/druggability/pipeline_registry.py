"""Project pipelines."""

from __future__ import annotations

from kedro.framework.project import find_pipelines
from kedro.pipeline import Pipeline
from druggability.pipelines.extract_pockets.pipeline import create_pipeline as create_extract_pockets_pipeline
from druggability.pipelines.pocket_ml.pipeline import create_pipeline as create_pocket_ml_pipeline


def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    pipelines = find_pipelines(raise_errors=True)
    # pipelines["__default__"] = sum(pipelines.values())
    pipelines["extract_pockets"] = create_extract_pockets_pipeline()
    pipelines["pocket_ml"] = create_pocket_ml_pipeline()
    return pipelines
