"""Kedro pipeline for scraping training data from RCSB PDB."""

from kedro.pipeline import Pipeline, node

from .nodes import run_scraper


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
        node(
            func=run_scraper,
            inputs=[
                "params:scraping.output_dir",
                "params:scraping.output_json",
                "params:scraping.years",
                "params:scraping.resolution_min",
                "params:scraping.resolution_max",
                "params:scraping.rfree_min",
                "params:scraping.rfree_max",
                "params:scraping.limit",
                "params:scraping.polite_delay",
                "params:scraping.debug",
            ],
            outputs="scraping_summary",
            name="run_scraper_node",
        ),
    ])
