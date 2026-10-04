from __future__ import annotations

import pandas as pd
import plotly.express as px

GROUP_COLORS = {"SjD": "#b23a48", "nonSjD": "#457b9d", "HV": "#6a994e"}


def cohort_bar(counts: pd.DataFrame):
    return px.bar(counts, x="disease_group", y="patients", color="disease_group",
                  color_discrete_map=GROUP_COLORS, text_auto=True, labels={"disease_group": "Diagnostic group"})


def endotype_heatmap(means: pd.DataFrame):
    return px.imshow(means, aspect="auto", color_continuous_scale="RdBu_r", color_continuous_midpoint=0,
                     labels={"x": "MSG module", "y": "Endotype", "color": "Mean z-score"})


def projection_scatter(points: pd.DataFrame, module: str):
    color = "disease_group" if "disease_group" in points else None
    fig = px.scatter(points, x="MSG eigengene", y="Blood projected score", color=color,
                     color_discrete_map=GROUP_COLORS, hover_data={"MSG eigengene": ":.3f", "Blood projected score": ":.3f",
                                                                  **({"disease_group": True} if color else {})},
                     title=f"MSG–Blood projection: {module}")
    # Patient identifiers are deliberately absent from traces/customdata.
    return fig
