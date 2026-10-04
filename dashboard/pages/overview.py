from pathlib import Path
import plotly.express as px
import streamlit as st
from components.figures import cohort_bar
from components.layout import intro, source_box
from data.loaders import load_tsv, optional_tsv
from data.transforms import cohort_counts


def render(root: Path) -> None:
    st.header("Overview / Paired cohort")
    intro("What is the composition and quality of our matched Blood–MSG cohort?",
          "Patient-level composition is calculated directly from the paired manifest produced by Stage 00.")
    manifest = load_tsv(root / "00_pairing/paired_manifest.tsv", ("patient_id", "disease_group"), ("patient_id",))
    counts = cohort_counts(manifest)
    st.metric("Matched patients", int(manifest.patient_id.nunique()))
    st.plotly_chart(cohort_bar(counts), use_container_width=True)
    cols = st.columns(2)
    for col, label, filename in zip(cols, ("Unmatched Blood", "Unmatched MSG"), ("unmatched_blood.tsv", "unmatched_msg.tsv")):
        optional = optional_tsv(root / "00_pairing" / filename)
        col.metric(label, "Not available" if optional is None else len(optional))
    pca = root / "01_differential_expression/paired_joint_PCA.png"
    if pca.is_file():
        st.subheader("Existing paired-cohort PCA")
        st.image(str(pca), use_container_width=True)
    else:
        st.info("The optional Stage 01 paired PCA is not available.")
    with st.expander("Pairing and QC"):
        st.write("Stage 00 validates tissue pairing and diagnostic-group concordance. Counts above are descriptive views; no scientific model is run here.")
    source_box("00_pairing/paired_manifest.tsv", "00_pairing/unmatched_blood.tsv", "00_pairing/unmatched_msg.tsv")
