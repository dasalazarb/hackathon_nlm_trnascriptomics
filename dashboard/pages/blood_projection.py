from pathlib import Path
import streamlit as st
from components.figures import projection_scatter
from components.layout import intro, source_box
from data.loaders import load_tsv, optional_tsv
from data.transforms import projection_points
from data.validators import common_projection_modules

METRICS = ("n_all", "rho_all", "fdr_all", "rho_SjD", "fdr_SjD", "rho_residualized_group", "fdr_residualized_group")


def render(root: Path) -> None:
    st.header("Blood Projection Explorer")
    intro("Which MSG-defined molecular programs remain concordant when projected to paired Blood samples?",
          "Select a module present in all three official Stage 02 outputs.")
    base = root / "02_wgcna_endotypes_projection"
    msg = load_tsv(base / "02_MSG_module_eigengenes.tsv", ("patient_id",), ("patient_id",))
    blood = load_tsv(base / "05_Blood_projected_module_scores.tsv", ("patient_id",), ("patient_id",))
    required = ("module",) + METRICS
    concordance = load_tsv(base / "05_MSG_Blood_projection_concordance.tsv", required, ("module",))
    modules = common_projection_modules(msg, blood, concordance)
    module = st.selectbox("Projected module", modules)
    manifest = optional_tsv(root / "00_pairing/paired_manifest.tsv", ("patient_id", "disease_group"))
    points = projection_points(msg, blood, module, manifest)
    st.plotly_chart(projection_scatter(points, module), use_container_width=True)
    row = concordance.loc[concordance.module.astype(str) == module].iloc[0]
    cols = st.columns(4)
    for index, metric in enumerate(METRICS):
        value = row[metric]
        cols[index % 4].metric(metric, "NA" if value != value else f"{value:.3g}")
    coverage_columns = [c for c in concordance.columns if "gene" in c.lower() or "coverage" in c.lower()]
    if coverage_columns:
        st.caption("Projection gene coverage (official output)")
        st.dataframe(row[["module", *coverage_columns]].to_frame().T, hide_index=True, use_container_width=True)
    order = st.selectbox("Sort module table by", ["fdr_all", "rho_all", "fdr_SjD"])
    display = concordance[["module", *METRICS, *coverage_columns]].sort_values(order, ascending=not order.startswith("rho"))
    st.dataframe(display.style.format(precision=3, na_rep="NA"), hide_index=True, use_container_width=True)
    st.info("MSG module weights were fixed before projection to Blood; Blood scores are not newly trained modules.")
    st.caption("rho_all may also reflect between-group differences. Within-SjD and group-residualized estimates provide complementary interpretations.")
    with st.expander("Methods / Limitations"):
        st.write("The scatter performs only a validated one-to-one patient join. Correlations, P values, FDR, and coverage are read from the official R output and are not recalculated here.")
    source_box("02_MSG_module_eigengenes.tsv", "05_Blood_projected_module_scores.tsv", "05_MSG_Blood_projection_concordance.tsv")
