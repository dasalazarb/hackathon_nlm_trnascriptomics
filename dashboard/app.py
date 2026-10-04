from pathlib import Path
import sys
import streamlit as st

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from data.paths import OutputDirectoryError, resolve_output_dir, run_label
from data.loaders import DataContractError
from pages import overview, endotypes, blood_projection

st.set_page_config(page_title="Blood–MSG Transcriptomics Explorer", page_icon="🧬", layout="wide")
st.title("Blood–MSG Transcriptomics Explorer")
st.caption("Molecular endotypes and cross-tissue projection in Sjögren's disease")
page = st.sidebar.radio("Navigate", ("Overview", "Endotypes", "Blood Projection"))
if st.sidebar.button("Reload outputs"):
    st.cache_data.clear()
    st.rerun()

try:
    output = resolve_output_dir()
    st.sidebar.caption(f"Run: {run_label(output)}")
    st.sidebar.warning("Unversioned output directory")
    {"Overview": overview.render, "Endotypes": endotypes.render, "Blood Projection": blood_projection.render}[page](output)
except (OutputDirectoryError, FileNotFoundError, DataContractError, ValueError) as exc:
    st.error(str(exc))
    st.info("Run the required R stages separately, then point OUTPUT_DIR at that single completed output tree.")
