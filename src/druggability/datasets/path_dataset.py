from kedro.io import AbstractDataset

class PathDataset(AbstractDataset):
    def __init__(self, filepath: str):
        self._filepath = filepath

    def _load(self) -> str:
        return self._filepath

    def _save(self, data) -> None:
        raise NotImplementedError("Saving not supported.")
    
    def _describe(self) -> dict:
        return {"filepath": self._filepath}