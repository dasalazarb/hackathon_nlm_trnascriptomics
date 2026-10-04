"""Strict, read-only TSV loaders for R pipeline products."""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import streamlit as st


class DataContractError(ValueError):
    pass


@st.cache_data(show_spinner=False)
def load_tsv(path: Path, required: tuple[str, ...] = (), unique: tuple[str, ...] = ()) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required pipeline output is missing: {path.name}")
    try:
        frame = pd.read_csv(path, sep="\t")
    except Exception as exc:
        raise DataContractError(f"Could not read {path.name} as a tab-delimited file: {exc}") from exc
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise DataContractError(f"{path.name} is missing required column(s): {', '.join(missing)}")
    for column in unique:
        if frame[column].isna().any() or frame[column].duplicated().any():
            raise DataContractError(f"{path.name}: {column} must be non-missing and unique")
    return frame


def optional_tsv(path: Path, required: tuple[str, ...] = ()) -> pd.DataFrame | None:
    return load_tsv(path, required) if path.is_file() else None
