import gzip
from typing import Any, Callable, Dict

def process_partitions(partitions: Dict[str, Callable[[], Any]]) -> Dict[str, str]:
    result = {}

    for partition_id, filepath in partitions.items():
        with gzip.open(filepath(), "rt", encoding="utf-8") as f:
            text = f.read()
        name = partition_id.split('-')[0]
        result[name] = text

    return result, []