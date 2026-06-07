"""Project pipelines."""

from __future__ import annotations

from kedro.framework.project import find_pipelines
from druggability.pipelines.protein_unzip import create_pipeline as create_protein_unzip_pipeline
from druggability.pipelines.protein_align import create_pipeline as create_protein_align_pipeline
from druggability.pipelines.compare_pockets import create_pipeline as create_compare_pockets_pipeline
from kedro.pipeline import Pipeline


def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    # pipelines = find_pipelines()  # Automatically find pipelines in the pipelines/ directory
    pipelines = {
        "unzip": create_protein_unzip_pipeline(),
        "align": create_protein_align_pipeline(),
        "compare": create_compare_pockets_pipeline(),
    } 
    pipelines["__default__"] = sum(pipelines.values())
    return pipelines
