from kedro.pipeline import Pipeline, node, pipeline
from .nodes import process_partitions


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
        node(
            func=process_partitions,
            inputs="raw_protein_cifs",
            outputs="intermediate_protein_cifs",
            name="unzip_cif_files_node",
        )
    ])