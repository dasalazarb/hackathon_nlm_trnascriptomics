import pandas as pd
import pytest
from components.figures import projection_scatter
from data.loaders import DataContractError
from data.transforms import projection_points
from data.validators import common_projection_modules


def fixtures():
    msg = pd.DataFrame({"patient_id": ["P_SYN_1", "P_SYN_2"], "MEblue": [0.2, -0.1]})
    blood = pd.DataFrame({"patient_id": ["P_SYN_2", "P_SYN_1"], "blue": [0.4, 0.3]})
    concordance = pd.DataFrame({"module": ["blue"]})
    return msg, blood, concordance


def test_me_prefix_alignment():
    msg, blood, concordance = fixtures()
    assert common_projection_modules(msg, blood, concordance) == ["blue"]


def test_scatter_pairing_and_hover_exclude_patient_id():
    msg, blood, _ = fixtures()
    groups = pd.DataFrame({"patient_id": ["P_SYN_1", "P_SYN_2"], "disease_group": ["SjD", "HV"]})
    points = projection_points(msg, blood, "blue", groups)
    assert points.loc[points.patient_id == "P_SYN_1", "Blood projected score"].item() == 0.3
    fig = projection_scatter(points, "blue")
    rendered = fig.to_json()
    assert "P_SYN_1" not in rendered and "P_SYN_2" not in rendered
    assert "patient_id" not in rendered


def test_many_to_many_pairing_is_rejected():
    msg, blood, _ = fixtures()
    blood = pd.concat([blood, blood.iloc[[0]]], ignore_index=True)
    with pytest.raises(DataContractError, match="one-to-one"):
        projection_points(msg, blood, "blue")
