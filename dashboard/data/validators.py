from __future__ import annotations

import pandas as pd
from .loaders import DataContractError


def msg_modules(frame: pd.DataFrame) -> dict[str, str]:
    modules = {column[2:]: column for column in frame.columns if column.startswith("ME") and len(column) > 2}
    if not modules:
        raise DataContractError("MSG eigengene output has no ME<module> columns")
    return modules


def common_projection_modules(msg: pd.DataFrame, blood: pd.DataFrame, concordance: pd.DataFrame) -> list[str]:
    mapping = msg_modules(msg)
    common = sorted(set(mapping) & (set(blood.columns) - {"patient_id"}) & set(concordance["module"].astype(str)))
    if not common:
        raise DataContractError("No aligned modules exist across MSG, Blood, and concordance outputs")
    return common


def one_to_one(left: pd.DataFrame, right: pd.DataFrame, context: str) -> pd.DataFrame:
    try:
        return left.merge(right, on="patient_id", how="inner", validate="one_to_one")
    except pd.errors.MergeError as exc:
        raise DataContractError(f"{context} requires one-to-one patient_id pairing: {exc}") from exc
