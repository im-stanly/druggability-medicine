"""A Random Forest that classifies surface points as pocket or not."""

import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.preprocessing import StandardScaler


@dataclass
class PocketClassifier:
    model: RandomForestClassifier | None = None
    scaler: StandardScaler | None = None
    feature_names: list[str] = field(default_factory=list)

    n_estimators: int = 200
    max_depth: int | None = 15
    class_weight: str = "balanced"
    random_state: int = 42

    def fit(self, X: np.ndarray, y: np.ndarray,
            feature_names: list[str] | None = None) -> "PocketClassifier":
        valid = (y == 0) | (y == 1)
        X, y = X[valid], y[valid]

        self.feature_names = feature_names or []
        self.scaler = StandardScaler()
        X = self.scaler.fit_transform(X)

        self.model = RandomForestClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            class_weight=self.class_weight,
            random_state=self.random_state,
            n_jobs=-1,
        )
        self.model.fit(X, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.model.predict(self.scaler.transform(X))

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        return self.model.predict_proba(self.scaler.transform(X))[:, 1]

    def score(self, X: np.ndarray, y: np.ndarray) -> dict:
        valid = (y == 0) | (y == 1)
        X, y = X[valid], y[valid]
        pred = self.predict(X)
        prob = self.predict_proba(X)
        return {
            "n": len(X),
            "accuracy": float(accuracy_score(y, pred)),
            "roc_auc": float(roc_auc_score(y, prob)),
            "f1": float(f1_score(y, pred, zero_division=0)),
        }

    def importances(self) -> list[dict]:
        if self.model is None:
            return []
        imp = self.model.feature_importances_
        names = self.feature_names or [f"f{i}" for i in range(len(imp))]
        return sorted(
            ({"feature": n, "importance": float(v)} for n, v in zip(names, imp)),
            key=lambda x: x["importance"], reverse=True,
        )

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump({
                "model": self.model, "scaler": self.scaler,
                "feature_names": self.feature_names,
            }, f)

    @classmethod
    def load(cls, path: str | Path) -> "PocketClassifier":
        with open(path, "rb") as f:
            data = pickle.load(f)
        c = cls()
        c.model = data["model"]
        c.scaler = data["scaler"]
        c.feature_names = data.get("feature_names", [])
        return c

    def __repr__(self) -> str:
        if self.model is None:
            return "PocketClassifier(untrained)"
        return f"PocketClassifier({self.n_estimators} trees, depth={self.max_depth})"
