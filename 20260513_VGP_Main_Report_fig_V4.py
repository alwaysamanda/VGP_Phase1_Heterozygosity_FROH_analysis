#### ==== Script for interactive figures for Gardiner et al. 2026 VGP ROH analysis ==== ####
## Date: 20260513 (May 13th, 2026)
## Author: Amanda Gardiner
## Version: 4
## GOAL: Create main figures as interactive HTML plots 
## NOTES: Want to see the species name, lineage, and IUCN status labels
#### ==== ==== ####


########## ========== LOAD IN NECESSARY PACKAGES ========== ##########
import sys
import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import date
from scipy import stats as scipy_stats
from itertools import combinations

########## ========== LOAD IN DATA AND VARIABLES ========== ##########
data              = pd.read_csv(sys.argv[1])
results_directory = sys.argv[2]
today             = date.today()

########## ========== FILTERING CRITERIA ========== ##########
## -- Filter out DD/NE species -- ##
data = data[data["IUCN"].notna()]
data = data[~data["IUCN"].isin(["DD/NE"])]

iucn_levels = ["CR", "EN", "VU", "NT", "LC"]
data["IUCN"] = pd.Categorical(data["IUCN"], categories=iucn_levels, ordered=True)

## -- Filter out invertebrates -- ##
data = data[~data["Extended_lineage"].isin(["Cyclostomes", "Other Deuterostomes"])]

## -- Normalize NROH by number of chromosomes -- ##
data["NROH_normalized"] = data["all_aln_NROH"] / data["num_auto_chromosomes"]

## -- Set categories for IDRisk -- ##
data['IDRisk_category'] = pd.cut(
    data['IDRisk_longplus'],
    bins=[-float('inf'), 0.05, 0.25, 0.5, float('inf')],
    labels=['Low risk', 'Moderate risk', 'High risk', 'Extreme risk'],
    right=True
)

## -- Filter for wild only species -- ##
data_wild = data[data["Combined_captivity_flag_primary"].isin(["no", "no likely"])].copy()

########## ========== COLOR PALETTES ========== ##########
iucn_colors = {
    "LC": "#60C659",
    "NT": "#CCE226",
    "VU": "#F9E814",
    "EN": "#FC7F3F",
    "CR": "#D81E05",
}

Marine_flag_colors = {
    "yes": "#1A7A8A", 
    "no":  "#C8C5BC"
}

microhabitat_colors = {
  "Ter": "#C7A96A", 
  "Fos": "#A0522D", 
  "Arb": "#4CA66B", 
  "Aqu": "#2196A6",
  "Aer": "#6BAED6"
}

idrisk_colors = {
  "Low risk":      "#FFD700",
  "Moderate risk": "#DAA520",
  "High risk":     "#B8600A",
  "Extreme risk":  "#8B1A1A"
}

########## ========== Functions ========== ##########
## Convert hex colors to rgba
def hex_to_rgba(hex_color, alpha=0.5):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

## Compute boxplot stats for each jitterplot, with mean instead of median
def compute_box_stats(df, value_col):
    def boxstats(group):
        x   = group[value_col].dropna()
        q1  = x.quantile(0.25)
        q3  = x.quantile(0.75)
        iqr = q3 - q1
        return pd.Series({
            "q1":     q1,
            "q3":     q3,
            "middle": float(x.median()),
            "mean":   float(x.mean()),
            "ymin":   float(max(x.min(), q1 - 1.5 * iqr)),
            "ymax":   float(min(x.max(), q3 + 1.5 * iqr)),
        })
    result = (
        df.groupby("IUCN", observed=True)
        .apply(boxstats, include_groups=False)
        .reset_index()
    )
    return result


def add_boxplot_and_jitter(
    fig, row, col,
    data_wild, box_stats,
    value_col, transform,
    iucn_levels, iucn_colors,
    y_pos, box_half,
    rng, hover_ylabel,
    show_legend,
):
    """
    Add manual box shapes + jitter scatter traces for one subplot panel.
    Returns the list of shapes to be attached to that subplot axis.
    """
    shapes  = []
    y_order = list(reversed(iucn_levels))

    for lvl in iucn_levels:

        if lvl not in box_stats["IUCN"].values:
            continue

        stats  = box_stats[box_stats["IUCN"] == lvl].iloc[0]
        grp    = data_wild[data_wild["IUCN"] == lvl][
                     [value_col, "species", "Extended_lineage"]
                 ].dropna(subset=[value_col])
        pts    = grp[value_col]
        fill_c = iucn_colors[lvl]

        if pts.empty:
            continue

        yc        = y_pos[lvl]
        y0        = yc - box_half
        y1        = yc + box_half
        whisker_w = box_half * 0.6

        ## Apply transform to stat values
        sq_q1     = transform(max(stats["q1"],     0))
        sq_q3     = transform(max(stats["q3"],     0))
        sq_middle = transform(max(stats["middle"], 0))
        sq_mean   = transform(max(stats["mean"],   0))
        sq_ymin   = transform(max(stats["ymin"],   0))
        sq_ymax   = transform(max(stats["ymax"],   0))

        ## Axis references differ per subplot row
        xref = "x" if row == 1 else f"x{row}"
        yref = "y" if row == 1 else f"y{row}"

        ## -- IQR filled rectangle -- ##
        shapes.append(dict(
            type="rect", xref=xref, yref=yref,
            x0=sq_q1, x1=sq_q3, y0=y0, y1=y1,
            fillcolor=hex_to_rgba(fill_c, 0.5),
            line=dict(color="black", width=2.2),
            layer="below",
        ))

        ## -- Median line -- ##
        fig.add_trace(go.Scatter(
            x=[sq_middle, sq_middle],
            y=[y0, y1],
            mode="lines",
            line=dict(color="black", width=2.5),
            showlegend=False,
            legendgroup=lvl,
            hoverinfo="none",
        ), row=row, col=col)

        ## -- Mean line -- ##
        fig.add_trace(go.Scatter(
            x=[sq_mean],
            y=[yc],
            mode="markers",
            marker=dict(
                symbol="diamond",
                size=12,
                color="white",
                line=dict(color="black", width=2.5),
            ),
            showlegend=False,
            legendgroup=lvl,
            hoverinfo="none",
        ), row=row, col=col)

        ## -- Lower whisker + cap -- ##
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymin, x1=sq_q1, y0=yc, y1=yc,
            line=dict(color="black", width=2.2),
        ))
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymin, x1=sq_ymin,
            y0=yc - whisker_w, y1=yc + whisker_w,
            line=dict(color="black", width=2.2),
        ))

        ## -- Upper whisker + cap -- ##
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_q3, x1=sq_ymax, y0=yc, y1=yc,
            line=dict(color="black", width=2.2),
        ))
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymax, x1=sq_ymax,
            y0=yc - whisker_w, y1=yc + whisker_w,
            line=dict(color="black", width=2.2),
        ))

        ## -- Legend proxy -- ##
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            name=lvl,
            marker=dict(
                symbol="square", size=12,
                color=hex_to_rgba(fill_c, 0.5),
                line=dict(color="black", width=1.5),
            ),
            showlegend=show_legend,
            legendgroup=lvl,
            hoverinfo="none",
        ), row=row, col=col)

        ## -- Jitter points -- ##
        jitter_offset = rng.uniform(-0.08, 0.08, size=len(pts))
        x_jittered    = transform(pts.values) + jitter_offset * 0.05
        y_jittered    = [yc] * len(pts)
        custom        = np.column_stack([
            pts.values,
            grp["species"].values,
            grp["Extended_lineage"].values,
        ])

        fig.add_trace(go.Scatter(
            x=x_jittered,
            y=y_jittered,
            mode="markers",
            name=lvl,
            marker=dict(
                symbol="circle", size=10,
                color=fill_c, opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
            showlegend=False,
            legendgroup=lvl,
            customdata=custom,
            hovertemplate=(
                f"<b>{lvl}</b><br>"
                f"Species: %{{customdata[1]}}<br>"
                f"Lineage: %{{customdata[2]}}<br>"
                f"{hover_ylabel}: %{{customdata[0]:.3f}}<extra></extra>"
            ),
        ), row=row, col=col)

    return shapes


########## ========== MAIN FIGURE 1 ========== ##########
#### ---- Shared axis settings ---- ####
iucn_levels  = ["CR", "EN", "VU", "NT", "LC"]
y_order      = list(reversed(iucn_levels))
y_pos        = {lvl: i for i, lvl in enumerate(y_order)}
y_tickvals   = [y_pos[lvl] for lvl in y_order]
box_half     = 0.25
rng          = np.random.default_rng(42)

#### ---- Panel A: Heterozygosity ---- ####
sqrt_transform = np.sqrt

het_breaks    = np.arange(0, 32, 2)
het_tickvals  = list(np.sqrt(het_breaks))
het_ticktext  = [str(v) for v in het_breaks]

box_stats_het = compute_box_stats(data_wild, "heterozygosity")

#### ---- Panel B: FROH (ROH >= 1% of chromosome) ---- ####
froh_col     = "1%+_aln_FROH_percent"

froh_breaks   = np.arange(0, 110, 10)
froh_tickvals = list(np.sqrt(froh_breaks))
froh_ticktext = [str(v) for v in froh_breaks]

box_stats_froh = compute_box_stats(data_wild, froh_col)

#### ---- Create subplot figure ---- ####
fig = make_subplots(
    rows=2, cols=1,
    shared_xaxes=False,
    vertical_spacing=0.12,
    subplot_titles=[
        "Heterozygosity per kb outside ROH \u2265100kb",
        "Inbreeding Coefficient (%) for ROH \u22651% of their chromosome",
    ],
)

#### ---- Add Panel A ---- ####
shapes_het = add_boxplot_and_jitter(
    fig=fig, row=1, col=1,
    data_wild=data_wild,
    box_stats=box_stats_het,
    value_col="heterozygosity",
    transform=sqrt_transform,
    iucn_levels=iucn_levels,
    iucn_colors=iucn_colors,
    y_pos=y_pos,
    box_half=box_half,
    rng=rng,
    hover_ylabel="Heterozygosity",
    show_legend=True,
)

#### ---- Add Panel B ---- ####
shapes_froh = add_boxplot_and_jitter(
    fig=fig, row=2, col=1,
    data_wild=data_wild,
    box_stats=box_stats_froh,
    value_col=froh_col,
    transform=sqrt_transform,
    iucn_levels=iucn_levels,
    iucn_colors=iucn_colors,
    y_pos=y_pos,
    box_half=box_half,
    rng=rng,
    hover_ylabel="FROH (%)",
    show_legend=False,
)

#### ---- Layout ---- ####
fig.update_layout(
    shapes=shapes_het + shapes_froh,
    height=900,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=80, r=40, t=80, b=60),
    legend=dict(
        title=dict(text="<b>IUCN Status</b>"),
        traceorder="normal",
        font=dict(size=12),
    ),
)

## -- Panel A axes -- ##
fig.update_xaxes(
    tickvals=het_tickvals,
    ticktext=het_ticktext,
    tickfont=dict(size=12),
    range=[0, np.sqrt(max(het_breaks))],
    zeroline=False,
    showgrid=True,
    gridcolor="#dddddd",
    row=1, col=1,
)
fig.update_yaxes(
    tickvals=y_tickvals,
    ticktext=y_order,
    tickfont=dict(size=13),
    range=[-0.6, len(iucn_levels) - 0.4],
    showgrid=False,
    row=1, col=1,
)

## -- Panel 2 axes -- ##
fig.update_xaxes(
    tickvals=froh_tickvals,
    ticktext=froh_ticktext,
    tickfont=dict(size=12),
    range=[0, np.sqrt(max(froh_breaks))],
    zeroline=False,
    showgrid=True,
    gridcolor="#dddddd",
    row=2, col=1,
)
fig.update_yaxes(
    tickvals=y_tickvals,
    ticktext=y_order,
    tickfont=dict(size=13),
    range=[-0.6, len(iucn_levels) - 0.4],
    showgrid=False,
    row=2, col=1,
)

#### ---- Export HTML ---- ####
output_dir  = os.path.join(results_directory, "Main_Figures")
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, f"{today}_Main_Figure_1_mean.html")

fig.write_html(output_path, full_html=True, include_plotlyjs="cdn")
print(f"Saved: {output_path}")
########## ========== END FIG 1 ========== ##########





########## ========== MAIN FIGURE 2 ========== ##########
## -- Filter: species with FROH_all ≥5% -- ##
subset_dat_fivepercent_wild = data_wild[data_wild["all_aln_FROH_percent"] >= 5].copy()

## -- Calculate equations -- ##
equation_labels = []
for lvl in iucn_levels:
    grp = subset_dat_fivepercent_wild[subset_dat_fivepercent_wild["IUCN"] == lvl].copy()
    grp = grp[["all_aln_FROH_percent", "NROH_normalized"]].dropna()
    if len(grp) < 2:
        continue

    sqrt_froh = np.sqrt(grp["all_aln_FROH_percent"].values)
    sqrt_nroh = np.sqrt(grp["NROH_normalized"].values)

    slope, intercept, r_value, _, _ = scipy_stats.linregress(sqrt_froh, sqrt_nroh)
    r2 = r_value ** 2

    sign_str = "+" if intercept >= 0 else "-"
    label    = f"y = {slope:.2f}x {sign_str} {abs(intercept):.2f}  (R² = {r2:.2f})"

    equation_labels.append({
        "IUCN":      lvl,
        "slope":     slope,
        "intercept": intercept,
        "r2":        r2,
        "label":     label,
    })

eq_df = pd.DataFrame(equation_labels)

## -- Plot -- ##
fig2 = go.Figure()

## x-axis / y-axis display breaks (original scale, rendered on sqrt axis)
x_breaks = np.arange(0, 110, 10)          # 0, 10, 20, … 100
y_breaks = [0, 25, 50, 100, 200, 300, 400, 500, 600]

x_tickvals = np.sqrt(x_breaks)
y_tickvals = np.sqrt(y_breaks)

## x range for regression line drawing
x_raw_max  = subset_dat_fivepercent_wild["all_aln_FROH_percent"].max()
x_line_raw = np.linspace(0, x_raw_max * 1.03, 300)
x_line_sq  = np.sqrt(x_line_raw)

for lvl in iucn_levels:
    grp = subset_dat_fivepercent_wild[subset_dat_fivepercent_wild["IUCN"] == lvl].copy()
    grp = grp[["all_aln_FROH_percent", "NROH_normalized", "species", "Extended_lineage"]].dropna(
        subset=["all_aln_FROH_percent", "NROH_normalized"]
    )
    if grp.empty:
        continue

    color = iucn_colors[lvl]

    ## ---- Regression line + CI band ---- ##
    eq_row = eq_df[eq_df["IUCN"] == lvl]
    if not eq_row.empty:
        slope     = eq_row["slope"].values[0]
        intercept = eq_row["intercept"].values[0]
        r2        = eq_row["r2"].values[0]
        label     = eq_row["label"].values[0]

        ## Fitted line in sqrt-sqrt space, then back to original for hover
        y_line_sq = slope * x_line_sq + intercept
        y_line_sq = np.maximum(y_line_sq, 0)           # clamp negatives

        ## Residual SE for the CI band
        grp_fit   = grp.copy()
        sqrt_froh = np.sqrt(grp_fit["all_aln_FROH_percent"].values)
        sqrt_nroh = np.sqrt(grp_fit["NROH_normalized"].values)
        n         = len(sqrt_froh)

        if n > 2:
            y_hat   = slope * sqrt_froh + intercept
            resid   = sqrt_nroh - y_hat
            se_res  = np.sqrt(np.sum(resid**2) / (n - 2))
            x_mean  = sqrt_froh.mean()
            ss_x    = np.sum((sqrt_froh - x_mean) ** 2)

            se_line = se_res * np.sqrt(
                1/n + (x_line_sq - x_mean)**2 / ss_x
            )
            ci_upper = np.maximum(y_line_sq + se_line, 0)
            ci_lower = np.maximum(y_line_sq - se_line, 0)
        else:
            ci_upper = y_line_sq
            ci_lower = y_line_sq

        ## CI band (filled)
        fig2.add_trace(go.Scatter(
            x=np.concatenate([x_line_sq, x_line_sq[::-1]]),
            y=np.concatenate([ci_upper,  ci_lower[::-1]]),
            fill="toself",
            fillcolor=hex_to_rgba(color, 0.15),
            line=dict(width=0),
            showlegend=False,
            legendgroup=lvl,
            hoverinfo="none",
        ))

        ## Regression line
        fig2.add_trace(go.Scatter(
            x=x_line_sq,
            y=y_line_sq,
            mode="lines",
            line=dict(color=color, width=2),
            showlegend=False,
            legendgroup=lvl,
            hoverinfo="none",
        ))

        ## Equation annotation — placed near the top of given IUCN category's data
        x_annot = np.sqrt(grp["all_aln_FROH_percent"].quantile(0.85))
        y_annot = np.sqrt(grp["NROH_normalized"].quantile(0.92))

        fig2.add_annotation(
            x=x_annot,
            y=y_annot,
            text=label,
            showarrow=False,
            font=dict(size=11, color=color),
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="rgba(0,0,0,0)",
            xanchor="left",
        )

    ## ---- Scatter points ---- ##
    custom = np.column_stack([
        grp["all_aln_FROH_percent"].values,
        grp["NROH_normalized"].values,
        grp["species"].values,
        grp["Extended_lineage"].values
    ])

    fig2.add_trace(go.Scatter(
        x=np.sqrt(grp["all_aln_FROH_percent"].values),
        y=np.sqrt(grp["NROH_normalized"].values),
        mode="markers",
        name=lvl,
        legendgroup=lvl,
        marker=dict(
            symbol="circle",
            size=12,
            color=color,
            opacity=0.8,
            line=dict(color="black", width=1.2),
        ),
        customdata=custom,
        hovertemplate=(
            f"<b>{lvl}</b><br>"
            "Species: %{customdata[2]}<br>"
            "Lineage: %{customdata[3]}<br>" 
            "FROH (%): %{customdata[0]:.2f}<br>"
            "NROH (norm.): %{customdata[1]:.2f}<extra></extra>"
        ),
        showlegend=True,
    ))

#### ---- Axes & layout ---- ####
fig2.update_layout(
    height=750,
    width=750,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=80, r=40, t=60, b=80),
    legend=dict(
        title=dict(text="<b>IUCN Status</b>"),
        font=dict(size=12),
    ),
    title=dict(
        text="F<sub>ROH</sub> (%) vs N<sub>ROH</sub> — species with F<sub>ROH</sub> ≥ 5%",
        font=dict(size=14),
        x=0.5,
    ),
)

fig2.update_xaxes(
    tickvals=list(x_tickvals),
    ticktext=[str(v) for v in x_breaks],
    tickfont=dict(size=12),
    title=dict(
        text="<b>F<sub>ROH</sub> (%) — ROH ≥ 100 kb (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, np.sqrt(x_raw_max * 1.05)],
    zeroline=False,
    showgrid=True,
    gridcolor="#dddddd",
)

fig2.update_yaxes(
    tickvals=list(y_tickvals),
    ticktext=[str(v) for v in y_breaks],
    tickfont=dict(size=12),
    title=dict(
        text="<b>N<sub>ROH</sub> (normalised, √ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, None],
    zeroline=False,
    showgrid=True,
    gridcolor="#dddddd",
)

#### ---- Export ---- ####
output_path_fig2 = os.path.join(output_dir, f"{today}_Main_Figure_2.html")
fig2.write_html(output_path_fig2, full_html=True, include_plotlyjs="cdn")
print(f"Saved: {output_path_fig2}")

########## ========== END FIG 2 ========== ##########





########## ========== MAIN FIGURE 3 ========== ##########
#### ---- Constants ---- ####
mu            = 1.25e-8
cutoff_old    = 3000
cutoff_recent = 1000

#### ---- Read MSMC files ---- ####
msmc_frames = []

for _, row in data_wild.iterrows():
    fpath   = row["MSMC_file_name"]
    species = row["species"]

    if not os.path.isfile(fpath):
        print(f"Skipping missing file: {fpath}")
        continue

    try:
        spec_dat = pd.read_table(fpath)
        spec_dat["right_time_boundary"] = pd.to_numeric(
            spec_dat["right_time_boundary"], errors="coerce"
        )

        msmc_frames.append(pd.DataFrame({
            "Generations_end":   spec_dat["left_time_boundary"]  / mu,
            "Generations_start": spec_dat["right_time_boundary"] / mu,
            "mid_gen":           (spec_dat["left_time_boundary"] + spec_dat["right_time_boundary"]) / (2 * mu),
            "Ne":                (1 / spec_dat["lambda"]) / (2 * mu),
            "Species":           species,
        }))
    except Exception as e:
        print(f"Error reading {fpath}: {e}")

data_read = pd.concat(msmc_frames, ignore_index=True)

## -- Attach IUCN + lineage from main data -- ##
meta_cols = data[["species", "IUCN", "Extended_lineage"]].drop_duplicates()
data_read = data_read.merge(meta_cols, left_on="Species", right_on="species", how="left")
data_read["IUCN"] = pd.Categorical(data_read["IUCN"], categories=iucn_levels, ordered=True)

#### ---- Build common log time grid ---- ####
valid_starts = data_read.loc[data_read["Generations_start"] > 0, "Generations_start"]
valid_ends   = data_read["Generations_end"]

time_grid = np.exp(
    np.linspace(
        np.log(valid_starts.min()),
        np.log(valid_ends.max()),
        200,
    )
)

#### ---- Interpolate each species onto the grid ---- ####
interp_frames = []

for (species, iucn, lineage), grp in data_read.groupby(
    ["Species", "IUCN", "Extended_lineage"], observed=True
):
    grp_sorted = grp.sort_values("mid_gen")
    ne_interp  = np.interp(
        time_grid,
        grp_sorted["mid_gen"].values,
        grp_sorted["Ne"].values,
    )
    interp_frames.append(pd.DataFrame({
        "Generations":      time_grid,
        "Ne":               ne_interp,
        "Species":          species,
        "IUCN":             iucn,
        "Extended_lineage": lineage,
    }))

data_interp = pd.concat(interp_frames, ignore_index=True)

#### ---- Summarise onto median / IQR per IUCN x time point ---- ####
data_summary = (
    data_interp
    .groupby(["IUCN", "Generations"], observed=True)["Ne"]
    .agg(
        median_Ne=lambda x: x.median(),
        lq_Ne    =lambda x: x.quantile(0.25),
        uq_Ne    =lambda x: x.quantile(0.75),
    )
    .reset_index()
)

## -- Save CSV (mirrors R write.csv) -- ##
csv_path = os.path.join(results_directory, f"{today}_wild_MSMC_grouped_stats.csv")
data_summary.to_csv(csv_path, index=False)
print(f"Saved summary CSV: {csv_path}")

#### ---- Jitterplot helper (reuses add_boxplot_and_jitter from Fig 1) ---- ####
## For Figure 3 panels B & C we need Ne on a log10 x-axis with IUCN on y.
## We reuse the existing helper but pass log10 as the transform.

log10_safe = lambda v: np.log10(np.maximum(v, 1))  # guard against 0 / negatives

## Log10 tick mapping for Ne axis
ne_breaks   = [100, 1_000, 10_000, 100_000, 1_000_000]
ne_tickvals = [log10_safe(v) for v in ne_breaks]
ne_ticktext = ["100", "1 k", "10 k", "100 k", "1 M"]

#### ---- Build Figure 3 (3-panel subplot) ---- ####
fig3 = make_subplots(
    rows=3, cols=1,
    subplot_titles=[
        "MSMC — Effective Population Size through time by IUCN status",
        f"Contemporary N<sub>e</sub>  (Generations end = 0)",
        f"Historical N<sub>e</sub>  ({cutoff_old:,} ≤ Generations ≤ start boundary)",
    ],
    vertical_spacing=0.10,
    row_heights=[0.50, 0.25, 0.25],
)

#### ---- Panel A: MSMC ribbon + line plot ---- ####
shown_in_legend = set()

for lvl in iucn_levels:
    sub = data_summary[data_summary["IUCN"] == lvl]
    if sub.empty:
        continue

    color       = iucn_colors[lvl]
    show_leg    = lvl not in shown_in_legend
    shown_in_legend.add(lvl)

    ## IQR ribbon
    fig3.add_trace(go.Scatter(
        x=np.concatenate([sub["Generations"], sub["Generations"][::-1]]),
        y=np.concatenate([sub["uq_Ne"],        sub["lq_Ne"][::-1]]),
        fill="toself",
        fillcolor=hex_to_rgba(color, 0.25),
        line=dict(width=0),
        showlegend=False,
        legendgroup=lvl,
        hoverinfo="none",
        name=lvl,
    ), row=1, col=1)

    ## Median line
    fig3.add_trace(go.Scatter(
        x=sub["Generations"],
        y=sub["median_Ne"],
        mode="lines",
        line=dict(color=color, width=2),
        name=lvl,
        legendgroup=lvl,
        showlegend=show_leg,
        hovertemplate=(
            f"<b>IUCN: {lvl}</b><br>"
            "Generations: %{x:,.0f}<br>"
            "Median N<sub>e</sub>: %{y:,.0f}<extra></extra>"
        ),
    ), row=1, col=1)

## Individual species lines (thin, semi-transparent) with hover labels
for (species, iucn, lineage), grp in data_interp.groupby(
    ["Species", "IUCN", "Extended_lineage"], observed=True
):
    color = iucn_colors.get(str(iucn), "#888888")
    fig3.add_trace(go.Scatter(
        x=grp["Generations"],
        y=grp["Ne"],
        mode="lines",
        line=dict(color=color, width=0.8),
        opacity=0.35,
        showlegend=False,
        legendgroup=str(iucn),
        name=species,
        hovertemplate=(
            f"<b>{species}</b><br>"
            f"IUCN: {iucn}<br>"
            f"Lineage: {lineage}<br>"
            "Generations: %{x:,.0f}<br>"
            "N<sub>e</sub>: %{y:,.0f}<extra></extra>"
        ),
    ), row=1, col=1)

## Vertical cutoff lines
for xval, label in [(cutoff_recent, f"{cutoff_recent:,} gen"),
                    (cutoff_old,    f"{cutoff_old:,} gen")]:
    fig3.add_vline(
        x=np.log10(xval),
        line=dict(color="black", width=1.5, dash="dash"),
        row=1, col=1,
    )
    fig3.add_annotation(
        x=np.log10(xval), y=1.02,
        xref="x", yref="paper",
        text=label,
        showarrow=False,
        font=dict(size=10),
        xanchor="left",
    )

#### ---- Panel B: Contemporary Ne (Generations_end == 0) ---- ####
modern_dat = data_read[data_read["Generations_end"] == 0].copy()
modern_dat["IUCN"] = pd.Categorical(modern_dat["IUCN"], categories=iucn_levels, ordered=True)

box_stats_modern = compute_box_stats(modern_dat, "Ne")

shapes_modern = add_boxplot_and_jitter(
    fig=fig3, row=2, col=1,
    data_wild=modern_dat,
    box_stats=box_stats_modern,
    value_col="Ne",
    transform=log10_safe,
    iucn_levels=iucn_levels,
    iucn_colors=iucn_colors,
    y_pos=y_pos,
    box_half=box_half,
    rng=rng,
    hover_ylabel="N<sub>e</sub>",
    show_legend=False,
)

#### ---- Panel C: Historical Ne (window spanning cutoff_old) ---- ####
old_dat = data_read[
    (data_read["Generations_end"] <= cutoff_old) &
    (data_read["Generations_start"] >= cutoff_old)
].copy()
old_dat["IUCN"] = pd.Categorical(old_dat["IUCN"], categories=iucn_levels, ordered=True)

box_stats_old = compute_box_stats(old_dat, "Ne")

shapes_old = add_boxplot_and_jitter(
    fig=fig3, row=3, col=1,
    data_wild=old_dat,
    box_stats=box_stats_old,
    value_col="Ne",
    transform=log10_safe,
    iucn_levels=iucn_levels,
    iucn_colors=iucn_colors,
    y_pos=y_pos,
    box_half=box_half,
    rng=rng,
    hover_ylabel="N<sub>e</sub>",
    show_legend=False,
)

#### ---- Axes & layout ---- ####
fig3.update_layout(
    shapes=shapes_modern + shapes_old,
    height=1200,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=90, r=40, t=80, b=60),
    legend=dict(
        title=dict(text="<b>IUCN Status</b>"),
        font=dict(size=12),
    ),
)

## Panel A — log10 x and y
fig3.update_xaxes(
    type="log",
    tickformat=".0e",
    title=dict(text="<b>Generations</b>", font=dict(size=13)),
    tickfont=dict(size=11),
    showgrid=True, gridcolor="#dddddd",
    zeroline=False,
    row=1, col=1,
)
fig3.update_yaxes(
    type="log",
    tickformat=",",
    title=dict(text="<b>Effective Population Size</b>", font=dict(size=13)),
    tickfont=dict(size=11),
    showgrid=True, gridcolor="#dddddd",
    zeroline=False,
    row=1, col=1,
)

## Panels B & C — log10 Ne on x, IUCN categories on y
for row_idx in [2, 3]:
    fig3.update_xaxes(
        tickvals=ne_tickvals,
        ticktext=ne_ticktext,
        tickfont=dict(size=12),
        title=dict(text="<b>N<sub>e</sub> (log<sub>10</sub> scale)</b>", font=dict(size=13)),
        showgrid=True, gridcolor="#dddddd",
        zeroline=False,
        row=row_idx, col=1,
    )
    fig3.update_yaxes(
        tickvals=y_tickvals,
        ticktext=y_order,
        tickfont=dict(size=13),
        range=[-0.6, len(iucn_levels) - 0.4],
        showgrid=False,
        row=row_idx, col=1,
    )

#### ---- Export ---- ####
output_path_fig3 = os.path.join(output_dir, f"{today}_Main_Figure_3.html")
fig3.write_html(output_path_fig3, full_html=True, include_plotlyjs="cdn")
print(f"Saved: {output_path_fig3}")

########## ========== END FIG 3 ========== ##########





########## ========== MAIN FIGURE 4 ========== ##########
## -- Prep data -- ##
data_wild = data_wild.rename(columns={"1%+_aln_FROH_percent": "longplus_froh"})

## -- Create columns for habitat count -- ##
def count_major_habitats(hab_str):
    if pd.isna(hab_str):
        return pd.NA
    parts   = [h.strip() for h in hab_str.split(";")]
    majors  = [h.split("_")[0] for h in parts if h]
    return len(set(majors))

def count_all_habitats(hab_str):
    if pd.isna(hab_str):
        return pd.NA
    return len([h for h in hab_str.split(";") if h.strip()])

wild_gen_hab_dat = data_wild[[
    "species", "IUCN", "Combined_captivity_flag_primary",
    "habitats", "longplus_froh", "all_aln_FROH_percent", 
    "heterozygosity", "Extended_lineage",  
]].copy()

wild_gen_hab_dat["major_habitat_count"] = wild_gen_hab_dat["habitats"].apply(count_major_habitats)
wild_gen_hab_dat["all_habitat_count"]   = wild_gen_hab_dat["habitats"].apply(count_all_habitats)

## -- Create the marine flag -- ##
def assign_marine_flag(hab_str):
    if pd.isna(hab_str):
        return pd.NA
    hab_str = hab_str.strip()
    return "yes" if re.search(r'(^|;\s*)(9|10|11|12)_', hab_str) else "no"

df_major_hab = wild_gen_hab_dat.copy()
df_major_hab["Marine_flag"] = df_major_hab["habitats"].apply(assign_marine_flag)
df_major_hab = df_major_hab[df_major_hab["Marine_flag"].notna()].copy()

## -- Group by microhabitat -- ##
microhabitat_cols = ["Fos", "Aer", "Aqu", "Ter", "Arb"]

data_subset_wild = data_wild[[
    "species", "IUCN", "Extended_lineage",
    "longplus_froh", "all_aln_FROH_percent", "heterozygosity",
] + microhabitat_cols].copy()

data_subset_long_wild = data_subset_wild.melt(
    id_vars=["species", "IUCN", "Extended_lineage",
             "longplus_froh", "all_aln_FROH_percent", "heterozygosity"],
    value_vars=microhabitat_cols,
    var_name="group",
    value_name="presence",
)
data_subset_long_wild = data_subset_long_wild[
    data_subset_long_wild["presence"] == 1
].copy()

microhabitat_order = ["Aer", "Arb", "Fos", "Ter", "Aqu"]
data_subset_long_wild["group"] = pd.Categorical(
    data_subset_long_wild["group"],
    categories=microhabitat_order,
    ordered=True,
)
data_subset_long_wild = data_subset_long_wild.sort_values("group")

## -- Statistical tests -- ##
def dunn_test_bh(df, value_col, group_col):
    """
    Pairwise Dunn test with Benjamini-Hochberg correction.
    Returns a DataFrame of significant pairs with annotations.
    """
    groups      = df[group_col].dropna().unique()
    pairs       = list(combinations(groups, 2))
    group_data  = [df.loc[df[group_col] == g, value_col].dropna().values for g in groups]

    # Kruskal-Wallis first
    kw_stat, kw_p = scipy_stats.kruskal(*group_data)
    print(f"  Kruskal-Wallis [{value_col} ~ {group_col}]: H={kw_stat:.3f}, p={kw_p:.4f}")

    # Pairwise Mann-Whitney U as Dunn approximation
    raw_results = []
    for (g1, g2) in pairs:
        x1 = df.loc[df[group_col] == g1, value_col].dropna().values
        x2 = df.loc[df[group_col] == g2, value_col].dropna().values
        if len(x1) < 2 or len(x2) < 2:
            continue
        _, p = scipy_stats.mannwhitneyu(x1, x2, alternative="two-sided")
        raw_results.append({"group1": g1, "group2": g2, "p_raw": p})

    if not raw_results:
        return pd.DataFrame()

    res_df   = pd.DataFrame(raw_results)
    n_tests  = len(res_df)
    # BH correction
    rank     = res_df["p_raw"].rank()
    res_df["p_adj"] = (res_df["p_raw"] * n_tests / rank).clip(upper=1.0)

    def annotate(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return None

    res_df["annotation"] = res_df["p_adj"].apply(annotate)
    sig = res_df[res_df["annotation"].notna()].copy()
    print(f"  Significant pairs: {len(sig)}")
    return sig

## -- Plot -- ##
def compute_box_stats_group(df, value_col, group_col):
    """Like compute_box_stats but for an arbitrary group column."""
    def boxstats(grp):
        x   = grp[value_col].dropna()
        q1  = x.quantile(0.25)
        q3  = x.quantile(0.75)
        iqr = q3 - q1
        return pd.Series({
            "q1":     q1,
            "q3":     q3,
            "middle": float(x.median()),
            "mean":   float(x.mean()),
            "ymin":   float(max(x.min(), q1 - 1.5 * iqr)),
            "ymax":   float(min(x.max(), q3 + 1.5 * iqr)),
        })
    return (
        df.groupby(group_col, observed=True)
        .apply(boxstats, include_groups=False)
        .reset_index()
        .rename(columns={group_col: "group"})
    )


def add_horizontal_boxplot_jitter(
    fig, row, col,
    data_df, box_stats,
    value_col, transform,
    group_order, group_colors,
    hover_ylabel,
    show_legend,
    rng,
    box_half=0.25,
    min_box_width=0.05,
    group_col="group",
):
    """
    Horizontal boxplot + jitter for an arbitrary categorical grouping.
    Categories on y-axis, transformed value on x-axis.
    """
    shapes   = []
    y_pos    = {g: i for i, g in enumerate(group_order)}

    for grp_label in group_order:
        row_stat = box_stats[box_stats["group"] == grp_label]
        if row_stat.empty:
            continue

        stats  = row_stat.iloc[0]
        pts_df = data_df[data_df[group_col] == grp_label][
            [value_col, "species", "IUCN", "Extended_lineage"]
        ].dropna(subset=[value_col])
        pts    = pts_df[value_col]

        if pts.empty:
            continue

        fill_c    = group_colors.get(grp_label, "#888888")
        yc        = y_pos[grp_label]
        y0        = yc - box_half
        y1        = yc + box_half
        whisker_w = box_half * 0.6

        sq_q1     = transform(max(stats["q1"],     0))
        sq_q3     = transform(max(stats["q3"],     0))
        sq_middle = transform(max(stats["middle"], 0))
        sq_mean   = transform(max(stats["mean"],   0))
        sq_ymin   = transform(max(stats["ymin"],   0))
        sq_ymax   = transform(max(stats["ymax"],   0))

        # Enforce minimum box width
        if (sq_q3 - sq_q1) < min_box_width:
            mid   = (sq_q1 + sq_q3) / 2
            sq_q1 = mid - min_box_width / 2
            sq_q3 = mid + min_box_width / 2

        xref = "x" if row == 1 else f"x{row}"
        yref = "y" if row == 1 else f"y{row}"

        ## IQR box
        shapes.append(dict(
            type="rect", xref=xref, yref=yref,
            x0=sq_q1, x1=sq_q3, y0=y0, y1=y1,
            fillcolor=hex_to_rgba(fill_c, 0.5),
            line=dict(color="black", width=2.2),
            layer="below",
        ))

        ## Median line
        fig.add_trace(go.Scatter(
            x=[sq_middle, sq_middle], y=[y0, y1],
            mode="lines",
            line=dict(color="black", width=2.5),
            showlegend=False, legendgroup=grp_label, hoverinfo="none",
        ), row=row, col=col)

        ## Mean diamond
        fig.add_trace(go.Scatter(
            x=[sq_mean], y=[yc],
            mode="markers",
            marker=dict(
                symbol="diamond", size=12,
                color="white", line=dict(color="black", width=2.5),
            ),
            showlegend=False, legendgroup=grp_label, hoverinfo="none",
        ), row=row, col=col)

        ## Lower whisker + cap
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymin, x1=sq_q1, y0=yc, y1=yc,
            line=dict(color="black", width=2.2),
        ))
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymin, x1=sq_ymin,
            y0=yc - whisker_w, y1=yc + whisker_w,
            line=dict(color="black", width=2.2),
        ))

        ## Upper whisker + cap
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_q3, x1=sq_ymax, y0=yc, y1=yc,
            line=dict(color="black", width=2.2),
        ))
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=sq_ymax, x1=sq_ymax,
            y0=yc - whisker_w, y1=yc + whisker_w,
            line=dict(color="black", width=2.2),
        ))

        ## Legend proxy
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            name=grp_label,
            marker=dict(
                symbol="square", size=12,
                color=hex_to_rgba(fill_c, 0.5),
                line=dict(color="black", width=1.5),
            ),
            showlegend=show_legend,
            legendgroup=grp_label, hoverinfo="none",
        ), row=row, col=col)

        ## Jitter points
        jitter_offset = rng.uniform(-0.08, 0.08, size=len(pts))
        x_jittered    = transform(pts.values) + jitter_offset * 0.05
        y_jittered    = [yc] * len(pts)
        custom        = np.column_stack([
            pts.values,
            pts_df["species"].values,
            pts_df["IUCN"].values,
            pts_df["Extended_lineage"].values,
        ])

        fig.add_trace(go.Scatter(
            x=x_jittered,
            y=y_jittered,
            mode="markers",
            name=grp_label,
            marker=dict(
                symbol="circle", size=10,
                color=fill_c, opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
            showlegend=False,
            legendgroup=grp_label,
            customdata=custom,
            hovertemplate=(
                f"<b>{grp_label}</b><br>"
                "Species: %{customdata[1]}<br>"
                "IUCN: %{customdata[2]}<br>"
                "Lineage: %{customdata[3]}<br>"
                f"{hover_ylabel}: %{{customdata[0]:.3f}}<extra></extra>"
            ),
        ), row=row, col=col)

    return shapes

## -- Add significance brackets -- ##
def add_significance_brackets(fig, sig_df, group_order, y_pos,
                               transform, x_max_transformed,
                               row, col, step=0.04):
    """
    Draw horizontal significance brackets above the plot area.
    Brackets extend past x_max_transformed, stepping outward per pair.
    """
    annotations = []
    shapes      = []
    xref        = "x" if row == 1 else f"x{row}"
    yref        = "y" if row == 1 else f"y{row}"

    for i, (_, pair_row) in enumerate(sig_df.iterrows()):
        g1    = pair_row["group1"]
        g2    = pair_row["group2"]
        annot = pair_row["annotation"]

        if g1 not in y_pos or g2 not in y_pos:
            continue

        y1    = y_pos[g1]
        y2    = y_pos[g2]
        x_tip = x_max_transformed * (1.04 + i * step)
        x_txt = x_max_transformed * (1.04 + i * step + 0.01)

        # Vertical bar at g1
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=x_tip, x1=x_tip, y0=y1, y1=y2,
            line=dict(color="black", width=1.5),
        ))
        # Tick at g1
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=x_tip * 0.995, x1=x_tip, y0=y1, y1=y1,
            line=dict(color="black", width=1.5),
        ))
        # Tick at g2
        shapes.append(dict(
            type="line", xref=xref, yref=yref,
            x0=x_tip * 0.995, x1=x_tip, y0=y2, y1=y2,
            line=dict(color="black", width=1.5),
        ))
        # Annotation at midpoint
        annotations.append(dict(
            x=x_txt,
            y=(y1 + y2) / 2,
            xref=xref, yref=yref,
            text=f"<b>{annot}</b>",
            showarrow=False,
            font=dict(size=13),
            xanchor="left",
        ))

    return shapes, annotations

## -- Run statistical tests -- ##
print("Marine flag ~ longplus_froh:")
sig_marine = dunn_test_bh(df_major_hab, "longplus_froh", "Marine_flag")

print("\nMicrohabitat ~ longplus_froh:")
sig_micro_froh = dunn_test_bh(data_subset_long_wild, "longplus_froh", "group")

print("\nMicrohabitat ~ heterozygosity:")
sig_micro_het = dunn_test_bh(data_subset_long_wild, "heterozygosity", "group")

## -- Axis ticks -- ##
sqrt_transform = np.sqrt

froh_breaks   = np.arange(0, 110, 10)
froh_tickvals = list(np.sqrt(froh_breaks))
froh_ticktext = [str(v) for v in froh_breaks]

het_breaks    = np.arange(0, 32, 2)
het_tickvals  = list(np.sqrt(het_breaks))
het_ticktext  = [str(v) for v in het_breaks]

marine_order      = ["yes", "no"]
microhabitat_order = ["Aer", "Arb", "Fos", "Ter", "Aqu"]

marine_y_pos = {g: i for i, g in enumerate(marine_order)}
micro_y_pos  = {g: i for i, g in enumerate(microhabitat_order)}

## -- Plot figure -- ##
fig4 = make_subplots(
    rows=3, cols=1,
    subplot_titles=[
        "F<sub>ROH</sub> (%) by Marine Habitat",
        "F<sub>ROH</sub> (%) by Microhabitat",
        "Heterozygosity by Microhabitat",
    ],
    vertical_spacing=0.10,
    row_heights=[0.25, 0.40, 0.40],
)

#### ---- Panel A: Marine flag ~ longplus_froh ---- ####
box_stats_marine = compute_box_stats_group(df_major_hab, "longplus_froh", "Marine_flag")

shapes_marine = add_horizontal_boxplot_jitter(
    fig=fig4, row=1, col=1,
    data_df=df_major_hab,
    box_stats=box_stats_marine,
    value_col="longplus_froh",
    transform=sqrt_transform,
    group_order=marine_order,
    group_colors=Marine_flag_colors,
    hover_ylabel="FROH (%)",
    show_legend=True,
    rng=rng,
    group_col="Marine_flag",
)

froh_marine_max = np.sqrt(df_major_hab["longplus_froh"].max())

if not sig_marine.empty:
    sig_marine_shapes, sig_marine_annots = add_significance_brackets(
        fig4, sig_marine,
        group_order=marine_order,
        y_pos=marine_y_pos,
        transform=sqrt_transform,
        x_max_transformed=froh_marine_max,
        row=1, col=1,
    )
else:
    sig_marine_shapes, sig_marine_annots = [], []

#### ---- Panel B: Microhabitat ~ longplus_froh ---- ####
box_stats_micro_froh = compute_box_stats_group(
    data_subset_long_wild, "longplus_froh", "group"
)

shapes_micro_froh = add_horizontal_boxplot_jitter(
    fig=fig4, row=2, col=1,
    data_df=data_subset_long_wild,
    box_stats=box_stats_micro_froh,
    value_col="longplus_froh",
    transform=sqrt_transform,
    group_order=microhabitat_order,
    group_colors=microhabitat_colors,
    hover_ylabel="FROH (%)",
    show_legend=False,
    rng=rng,
)

froh_micro_max = np.sqrt(data_subset_long_wild["longplus_froh"].max())

if not sig_micro_froh.empty:
    sig_micro_froh_shapes, sig_micro_froh_annots = add_significance_brackets(
        fig4, sig_micro_froh,
        group_order=microhabitat_order,
        y_pos=micro_y_pos,
        transform=sqrt_transform,
        x_max_transformed=froh_micro_max,
        row=2, col=1,
    )
else:
    sig_micro_froh_shapes, sig_micro_froh_annots = [], []

#### ---- Panel C: Microhabitat ~ heterozygosity ---- ####
box_stats_micro_het = compute_box_stats_group(
    data_subset_long_wild, "heterozygosity", "group"
)

shapes_micro_het = add_horizontal_boxplot_jitter(
    fig=fig4, row=3, col=1,
    data_df=data_subset_long_wild,
    box_stats=box_stats_micro_het,
    value_col="heterozygosity",
    transform=sqrt_transform,
    group_order=microhabitat_order,
    group_colors=microhabitat_colors,
    hover_ylabel="Heterozygosity",
    show_legend=False,
    rng=rng,
)

het_micro_max = np.sqrt(data_subset_long_wild["heterozygosity"].max())

if not sig_micro_het.empty:
    sig_micro_het_shapes, sig_micro_het_annots = add_significance_brackets(
        fig4, sig_micro_het,
        group_order=microhabitat_order,
        y_pos=micro_y_pos,
        transform=sqrt_transform,
        x_max_transformed=het_micro_max,
        row=3, col=1,
    )
else:
    sig_micro_het_shapes, sig_micro_het_annots = [], []

## -- Layout -- ##
all_shapes = (
    shapes_marine      + sig_marine_shapes +
    shapes_micro_froh  + sig_micro_froh_shapes +
    shapes_micro_het   + sig_micro_het_shapes
)
all_annotations = (
    sig_marine_annots +
    sig_micro_froh_annots +
    sig_micro_het_annots
)

fig4.update_layout(
    shapes=all_shapes,
    annotations=fig4.layout.annotations + tuple(all_annotations),
    height=1100,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=80, r=120, t=80, b=60),
    legend=dict(
        title=dict(text="<b>Group</b>"),
        orientation="h",
        yanchor="bottom", y=-0.08,
        xanchor="center", x=0.5,
        font=dict(size=12),
    ),
)

## -- Panel A axes -- ##
fig4.update_xaxes(
    tickvals=froh_tickvals,
    ticktext=froh_ticktext,
    tickfont=dict(size=12),
    title=dict(
        text="<b>F<sub>ROH</sub> (%) — ROH ≥1% chromosome (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, froh_marine_max * 1.25],
    showgrid=True, gridcolor="#dddddd", zeroline=False,
    row=1, col=1,
)
fig4.update_yaxes(
    tickvals=[marine_y_pos[g] for g in marine_order],
    ticktext=marine_order,
    tickfont=dict(size=13),
    range=[-0.6, len(marine_order) - 0.4],
    showgrid=False,
    row=1, col=1,
)

## -- Panel B axes -- ##
fig4.update_xaxes(
    tickvals=froh_tickvals,
    ticktext=froh_ticktext,
    tickfont=dict(size=12),
    title=dict(
        text="<b>F<sub>ROH</sub> (%) — ROH ≥1% chromosome (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, froh_micro_max * 1.30],
    showgrid=True, gridcolor="#dddddd", zeroline=False,
    row=2, col=1,
)
fig4.update_yaxes(
    tickvals=[micro_y_pos[g] for g in microhabitat_order],
    ticktext=microhabitat_order,
    tickfont=dict(size=13),
    range=[-0.6, len(microhabitat_order) - 0.4],
    showgrid=False,
    row=2, col=1,
)

## -- Panel C axes -- ##
fig4.update_xaxes(
    tickvals=het_tickvals,
    ticktext=het_ticktext,
    tickfont=dict(size=12),
    title=dict(
        text="<b>Heterozygosity per kb outside ROH (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, het_micro_max * 1.30],
    showgrid=True, gridcolor="#dddddd", zeroline=False,
    row=3, col=1,
)
fig4.update_yaxes(
    tickvals=[micro_y_pos[g] for g in microhabitat_order],
    ticktext=microhabitat_order,
    tickfont=dict(size=13),
    range=[-0.6, len(microhabitat_order) - 0.4],
    showgrid=False,
    row=3, col=1,
)

## -- Export -- ##
output_path_fig4 = os.path.join(output_dir, f"{today}_Main_Figure_4.html")
fig4.write_html(output_path_fig4, full_html=True, include_plotlyjs="cdn")
print(f"Saved: {output_path_fig4}")
########## ========== END FIG 4 ========== ##########





########## ========== MAIN FIGURE 5 ========== ##########
## -- Set shapes for IUCN categories -- ##
iucn_shapes = {
    "CR": "circle",
    "EN": "diamond",
    "VU": "triangle-up",
    "NT": "cross",
    "LC": "x-dot",
}

## -- Define IDRisk clines -- ##
cline_values = [0.05, 0.25, 0.5]
cline_colors = {
    0.05: "#DAA520",
    0.25: "#B8600A",
    0.50: "#8B1A1A",
}
cline_labels = {
    0.05: "IDRisk = 0.05",
    0.25: "IDRisk = 0.25",
    0.50: "IDRisk = 0.50",
}

froh_prop = data_wild["longplus_froh"] / 100
x_min     = froh_prop.min(skipna=True)
x_max     = froh_prop.max(skipna=True)

cline_frames = []
for k in cline_values:
    x_seq = np.linspace(x_min, x_max, 200)
    y_seq = k / x_seq          # y = k / x  →  IDRisk = FROH * het = k
    cline_frames.append(pd.DataFrame({"x": x_seq, "y": y_seq, "k": k}))
cline_dat = pd.concat(cline_frames, ignore_index=True)

## -- Label IDRisk regions defined by clines -- ##
label_x = x_max * 0.85

region_labels = pd.DataFrame({
    "x": [label_x] * 4,
    "y": [
        0.5  * cline_values[0] / label_x,
        np.sqrt(cline_values[0] / label_x * cline_values[1] / label_x),
        np.sqrt(cline_values[1] / label_x * cline_values[2] / label_x),
        1.5  * cline_values[2] / label_x,
    ],
    "label": ["Low risk", "Moderate risk", "High risk", "Extreme risk"],
    "color": ["#FFD700", "#DAA520", "#B8600A", "#8B1A1A"],
})

## -- Name outliers in data with very high IDRisk -- ##
data_wild["outlier_IDRisk_longplus"] = np.where(
    data_wild["IDRisk_longplus"] >= 0.5,
    data_wild["species"],
    None,
)

## -- Axis tick definition -- ##
froh_prop_breaks  = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
froh_tickvals     = list(np.sqrt(froh_prop_breaks))
froh_ticktext     = [f"{int(v*100)}%" for v in froh_prop_breaks]

het_breaks        = np.arange(0, 36, 5)
het_tickvals      = list(np.sqrt(het_breaks))
het_ticktext      = [str(v) for v in het_breaks]

## -- Build Fig 5 -- ##
fig5 = go.Figure()

## ---- Cline lines ---- ##
for k in cline_values:
    sub = cline_dat[cline_dat["k"] == k]
    fig5.add_trace(go.Scatter(
        x=np.sqrt(sub["x"]),
        y=np.sqrt(sub["y"]),
        mode="lines",
        line=dict(color=cline_colors[k], width=1.5, dash="solid"),
        name=cline_labels[k],
        legendgroup=f"cline_{k}",
        showlegend=True,
        hovertemplate=(
            f"IDRisk = {k}<br>"
            "FROH: %{customdata[0]:.3f}<br>"
            "Het: %{customdata[1]:.3f}<extra></extra>"
        ),
        customdata=np.column_stack([sub["x"], sub["y"]]),
    ))

## ---- Region band labels ---- ##
for _, lrow in region_labels.iterrows():
    fig5.add_annotation(
        x=np.sqrt(lrow["x"]),
        y=np.sqrt(lrow["y"]),
        text=f"<b>{lrow['label']}</b>",
        showarrow=False,
        font=dict(size=11, color=lrow["color"]),
        xanchor="right",
        bgcolor="rgba(255,255,255,0.7)",
        bordercolor="rgba(0,0,0,0)",
    )

## ---- Scatter points per IUCN level ---- ##
for lvl in iucn_levels:
    sub = data_wild[data_wild["IUCN"] == lvl].copy()
    sub = sub.dropna(subset=["longplus_froh", "heterozygosity", "IDRisk_longplus"])
    if sub.empty:
        continue

    marker_colors = sub["IDRisk_category"].map(idrisk_colors).fillna("#888888")

    custom = np.column_stack([
        sub["species"].values,
        sub["IUCN"].values,
        sub["Extended_lineage"].values,
        (sub["longplus_froh"] / 100).values,
        sub["heterozygosity"].values,
        sub["IDRisk_longplus"].values,
        sub["IDRisk_category"].astype(str).values,
    ])

    fig5.add_trace(go.Scatter(
        x=np.sqrt(sub["longplus_froh"] / 100),
        y=np.sqrt(sub["heterozygosity"]),
        mode="markers",
        name=lvl,
        legendgroup=lvl,
        marker=dict(
            symbol=iucn_shapes[lvl],
            size=10,
            color=marker_colors,
            line=dict(color="black", width=1.0),
        ),
        customdata=custom,
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "IUCN: %{customdata[1]}<br>"
            "Lineage: %{customdata[2]}<br>"
            "FROH (%): %{customdata[3]:.3f}<br>"
            "Heterozygosity: %{customdata[4]:.3f}<br>"
            "IDRisk: %{customdata[5]:.3f}<br>"
            "Risk category: %{customdata[6]}<extra></extra>"
        ),
        showlegend=True,
    ))

## ---- Outlier text labels ---- ##
outliers = data_wild[data_wild["outlier_IDRisk_longplus"].notna()].copy()
outliers = outliers.dropna(subset=["longplus_froh", "heterozygosity"])

if not outliers.empty:
    fig5.add_trace(go.Scatter(
        x=np.sqrt(outliers["longplus_froh"] / 100),
        y=np.sqrt(outliers["heterozygosity"]),
        mode="text",
        text=outliers["species"],
        textposition="top center",
        textfont=dict(size=9, color="black"),
        showlegend=False,
        hoverinfo="skip",
    ))

## ---- IDRisk color legend proxies ---- ##
for cat, col in idrisk_colors.items():
    fig5.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="markers",
        name=cat,
        legendgroup=cat,
        marker=dict(
            symbol="circle", size=10,
            color=col,
            line=dict(color="black", width=1),
        ),
        showlegend=True,
    ))

#### ---- Layout ---- ####
fig5.update_layout(
    height=750,
    width=800,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=80, r=40, t=60, b=80),
    legend=dict(
        title=dict(text="<b>IUCN / IDRisk</b>"),
        font=dict(size=11),
        tracegroupgap=8,
    ),
    title=dict(
        text="ID<sub>risk</sub> — F<sub>ROH</sub> vs Heterozygosity (wild species)",
        font=dict(size=14),
        x=0.5,
    ),
)

fig5.update_xaxes(
    tickvals=froh_tickvals,
    ticktext=froh_ticktext,
    tickfont=dict(size=12),
    title=dict(
        text="<b>F<sub>ROH</sub> for ROH ≥1% of chromosome (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, np.sqrt(x_max) * 1.05],
    showgrid=True,
    gridcolor="#dddddd",
    zeroline=False,
    showline=True,
    linecolor="black",
    mirror=True,
)

fig5.update_yaxes(
    tickvals=het_tickvals,
    ticktext=het_ticktext,
    tickfont=dict(size=12),
    title=dict(
        text="<b>Het/kb in non-ROH regions (√ scale)</b>",
        font=dict(size=13),
    ),
    range=[0, np.sqrt(35)],
    showgrid=True,
    gridcolor="#dddddd",
    zeroline=False,
    showline=True,
    linecolor="black",
    mirror=True,
)

#### ---- Export ---- ####
output_path_fig5 = os.path.join(output_dir, f"{today}_Main_Figure_5.html")
fig5.write_html(output_path_fig5, full_html=True, include_plotlyjs="cdn")
print(f"Saved: {output_path_fig5}")

########## ========== END FIG 5 ========== ##########