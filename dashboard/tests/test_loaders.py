from pathlib import Path
import pytest
from data.loaders import DataContractError, load_tsv


def test_loads_tab_delimited_synthetic_data(tmp_path: Path):
    path = tmp_path / "synthetic.tsv"
    path.write_text("patient_id\tdisease_group\nP_SYN_1\tSjD\n")
    assert load_tsv(path, ("patient_id",)).iloc[0].disease_group == "SjD"


def test_missing_required_column_is_readable(tmp_path: Path):
    path = tmp_path / "synthetic.tsv"
    path.write_text("wrong\nvalue\n")
    with pytest.raises(DataContractError, match="missing required column.*patient_id"):
        load_tsv(path, ("patient_id",))


def test_unique_patient_contract(tmp_path: Path):
    path = tmp_path / "synthetic.tsv"
    path.write_text("patient_id\nP_SYN_1\nP_SYN_1\n")
    with pytest.raises(DataContractError, match="must be non-missing and unique"):
        load_tsv(path, ("patient_id",), ("patient_id",))
