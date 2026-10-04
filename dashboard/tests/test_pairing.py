import pandas as pd
from data.transforms import cohort_counts, endotype_means


def test_cohort_counts_are_dynamic():
    synthetic = pd.DataFrame({"patient_id": ["P1", "P2", "P3"], "disease_group": ["SjD", "SjD", "HV"]})
    got = dict(cohort_counts(synthetic).set_index("disease_group").patients)
    assert got == {"HV": 1, "SjD": 2}


def test_missing_endotype_is_not_biological_category():
    assignments = pd.DataFrame({"patient_id": ["P1", "P2", "P3"], "endotype": ["E1", None, "E2"]})
    eig = pd.DataFrame({"patient_id": ["P1", "P2", "P3"], "MEblue": [1.0, 99.0, -1.0]})
    means = endotype_means(assignments, eig)
    assert list(means.index) == ["E1", "E2"]
    assert 99.0 not in means.blue.tolist()
