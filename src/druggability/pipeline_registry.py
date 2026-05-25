"""Project pipelines."""

from __future__ import annotations

from kedro.framework.project import find_pipelines
from druggability.pipelines.protein_unzip import create_pipeline as create_protein_unzip_pipeline
from druggability.pipelines.protein_compare import create_pipeline as create_protein_compare_pipeline
from kedro.pipeline import Pipeline


def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    pipelines = {}
    pipelines["unzip"] = create_protein_unzip_pipeline()
    pipelines["compare"] = create_protein_compare_pipeline()
    return pipelines
