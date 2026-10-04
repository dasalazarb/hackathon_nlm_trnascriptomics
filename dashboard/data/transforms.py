from __future__ import annotations

import pandas as pd
from .validators import msg_modules, one_to_one


def cohort_counts(manifest: pd.DataFrame) -> pd.DataFrame:
    patients = manifest[["patient_id", "disease_group"]].drop_duplicates()
    if patients["patient_id"].duplicated().any():
        raise ValueError("A patient has conflicting disease groups")
    return patients.groupby("disease_group", dropna=False).size().rename("patients").reset_index()


def endotype_means(assignments: pd.DataFrame, eigengenes: pd.DataFrame) -> pd.DataFrame:
    assigned = assignments.dropna(subset=["endotype"])
    joined = one_to_one(assigned[["patient_id", "endotype"]], eigengenes, "Endotype/eigengene join")
    columns = list(msg_modules(eigengenes).values())
    return joined.groupby("endotype", sort=True)[columns].mean().rename(columns=lambda c: c[2:])


def projection_points(msg: pd.DataFrame, blood: pd.DataFrame, module: str, groups: pd.DataFrame | None = None) -> pd.DataFrame:
    msg_column = f"ME{module}"
    if msg_column not in msg or module not in blood:
        raise ValueError(f"Module alignment failed: {msg_column} must map to {module}")
    points = one_to_one(msg[["patient_id", msg_column]], blood[["patient_id", module]], "MSG/Blood join")
    points = points.rename(columns={msg_column: "MSG eigengene", module: "Blood projected score"})
    if groups is not None:
        points = one_to_one(points, groups[["patient_id", "disease_group"]], "Diagnostic group join")
    return points
