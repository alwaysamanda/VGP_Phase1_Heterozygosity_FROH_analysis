#### ==== Script to calculate Linear Regression Model to predict ROH ==== ####
## Date: 20260527 (May 27th, 2026)
## Author: Amanda Gardiner
## Version 3 (V2 is 20260327_ROH_Linear_Regression_V2.py)
## GOAL: Create a linear regression model to see what predicts ROH for a species
## UPDATES: 20260618 -- Adding in linear model for heterozygosity and svg figures to save for main text
#### ====  ==== ####



########## ========== Load in necessary packages ========== ##########
from datetime import date
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import re as _re
import sys
from sys import stdin
import io
import os
from ete3 import Tree
from statsmodels.regression.linear_model import GLS
from scipy.special import logit, expit
import plotly.graph_objects as go
import plotly.express as px


########## ========== Load in data and variables ========== ##########
args = sys.argv
today = date.today()
data = pd.read_csv(args[1])
tree = Tree(args[2])
output_dir = args[3]


########## ========== Clean and prep data ========== ##########
## -- Only wild species -- ##
species_order = data.sort_values('clade')['species'].tolist()
data = data[~data['IUCN'].str.strip().str.upper().isin(['DD/NE', 'EW'])]
data_wild = data[data['Combined_captivity_flag_primary'].isin(['no', 'no likely'])]

## -- Want species with at least 1% FROH -- ##
data_wild = data_wild[
    (data_wild['all_aln_FROH_percent'] >= 1)
]

## -- Remove Na values -- ##
data_clean = data_wild.dropna(
    subset=[
        "Extended_lineage", 
        "IUCN", 
        "heterozygosity", 
        "all_aln_FROH_percent", 
        "habitats", 
        "Aqu", 
        "Fos", 
        "Ter", 
        "Arb", 
        "Aer"
        ], 
    how="any"
)

## -- Remove species with 0 in het or FROH -- ##
data_clean = data_clean[
    (data_clean['heterozygosity'] != 0) & 
    (data_clean['all_aln_FROH_percent'] != 0)
]

## -- Extract major habitat code and encode dummy variables for it -- ##
def extract_habitats(habitat_str):
    habitats = str(habitat_str).split(',')
    major_habitats = []
    for h in habitats:
        h = h.strip()
        if h:
            major = h.split('_')[0]
            major = major.split(';')[0].strip()
            if major:
                major_habitats.append(major)
    return list(set(major_habitats)) 

data_clean = data_clean.copy() 
data_clean['habitat_major'] = data_clean['habitats'].apply(extract_habitats)

mlb = MultiLabelBinarizer()
habitat_dummies = pd.DataFrame(
    mlb.fit_transform(data_clean['habitat_major']),
    columns=[f'hab_{c}' for c in mlb.classes_],
    index=data_clean.index
).astype(int)

data_clean = pd.concat([data_clean, habitat_dummies], axis=1)

## -- Extract IUCN status and encode as dummy variables -- ##
iucn_order = {'LC': 0, 'NT': 1, 'VU': 2, 'EN': 3, 'CR': 4,}
data_clean['IUCN_encoded'] = data_clean['IUCN'].map(iucn_order)


########## ========== Get covariance and run FROH model ========== ##########
## -- Match species in data to tips in tree -- ##
def build_phylo_cov_matrix(tree, accession_list):
    n = len(accession_list)
    cov = np.zeros((n, n))
    root = tree.get_tree_root()
    for i, acc1 in enumerate(accession_list):
        for j, acc2 in enumerate(accession_list):
            if i == j:
                node = tree.search_nodes(name=acc1)
                if not node:
                    raise ValueError(f"Accession '{acc1}' not found in tree")
                cov[i, j] = root.get_distance(node[0])
            else:
                ancestor = tree.get_common_ancestor(acc1, acc2)
                cov[i, j] = root.get_distance(ancestor)
    return pd.DataFrame(cov, index=accession_list, columns=accession_list)

tree_species = {leaf.name for leaf in tree.get_leaves()}
data_clean = data_clean[data_clean['primary_accession_num'].isin(tree_species)]
species_list = data_clean['species'].tolist()

## -- Build matrix and align to data -- ##
accession_list = data_clean['primary_accession_num'].tolist()
cov_matrix = build_phylo_cov_matrix(tree, accession_list)
cov_matrix.index = species_list
cov_matrix.columns = species_list
cov_array = cov_matrix.loc[species_list, species_list].values

## -- Fit PGLS -- ##
hab_cols = [c for c in data_clean.columns if c.startswith('hab_')]
fixed_effects = ['IUCN_encoded', 'heterozygosity', 'Aqu', 'Fos', 'Ter', 'Arb', 'Aer'] + hab_cols
X = sm.add_constant(data_clean[fixed_effects].astype(float))
y = data_clean['all_aln_FROH_percent'].astype(float)
y_logit = pd.Series(
    logit(y.values / 100).flatten(),
    index=y.index
)

model = GLS(y_logit, X, sigma=cov_array)
result = model.fit()
print(result.summary())


########## ========== Plot FROH results ========== ##########
#### ---- IUCN color and shape mappings ---- ####
iucn_colors = {
    "LC": "#60C659",
    "NT": "#CCE226",
    "VU": "#F9E814",
    "EN": "#FC7F3F",
    "CR": "#D81E05"
}

iucn_symbols = {
    "LC": "triangle-down",
    "NT": "triangle-up",
    "VU": "diamond",
    "EN": "square",
    "CR": "circle"
}

iucn_order = ["LC", "NT", "VU", "EN", "CR"]

FONT = "Arial"

#### ---- Convert habitat number to name ---- ####
habitat_name_map = {
    1:  "Forest",
    2:  "Savanna",
    3:  "Shrubland",
    4:  "Grassland",
    5:  "Wetlands (inland)",
    6:  "Rocky Areas",
    7:  "Caves (non-aquatic)",
    8:  "Desert",
    9:  "Marine Neritic",
    10: "Marine Oceanic",
    11: "Marine Deep Ocean Floor",
    12: "Marine Intertidal",
    13: "Marine Coastal/Supratidal",
    14: "Artificial Terrestrial",
    15: "Artificial Aquatic",
    16: "Introduced Vegetation",
    17: "Other",
}

def rename_coeff_label(label):
    m = _re.search(r'(\d+)$', str(label))  
    if m:
        num = int(m.group(1))
        if num in habitat_name_map:
            return habitat_name_map[num]
    return label


#### ---- Plot actual vs. predicted for FROH model ---- ####
y_froh_pred   = expit(result.fittedvalues) * 100
y_froh_actual = y.values

min_val = min(y_froh_actual.min(), y_froh_pred.min())
max_val = max(y_froh_actual.max(), y_froh_pred.max())

fig1 = go.Figure()

## -- Plot 1:1 reference line -- ##
fig1.add_trace(go.Scatter(
    x=[min_val, max_val],
    y=[min_val, max_val],
    mode='lines',
    line=dict(color='#888888', dash='dash', width=1.5),
    name='1 : 1 line',
    hoverinfo='skip'
))

## -- Add IUCN traces -- ##
for status in iucn_order:
    mask = data_clean['IUCN'] == status
    if not mask.any():
        continue

    subset     = data_clean[mask]
    actual_sub = y_froh_actual[mask]
    pred_sub   = y_froh_pred[mask]

    hover_sub = [
        f"<b>{row['species']}</b><br>"
        f"IUCN: {row['IUCN']}<br>"
        f"Lineage: {row['Extended_lineage']}<br>"
        f"FROH: {row['all_aln_FROH_percent']:.2f}%"
        for _, row in subset.iterrows()
    ]

    fig1.add_trace(go.Scatter(
        x=actual_sub,
        y=pred_sub,
        mode='markers',
        marker=dict(
            color=iucn_colors[status],
            symbol=iucn_symbols[status],
            size=10,
            opacity=0.9,
            line=dict(color='rgba(0,0,0,0.45)', width=0.8)
        ),
        text=hover_sub,
        hovertemplate="%{text}<br><b>Actual:</b> %{x:.2f}%<br><b>Predicted:</b> %{y:.2f}%<extra></extra>",
        name=status
    ))

## -- R² annotation -- ##
fig1.add_annotation(
    xref='paper', yref='paper',
    x=0.05, y=0.95,
    text=f"<b>R² = {result.rsquared:.4f}</b>",
    showarrow=False,
    bgcolor='white',
    bordercolor='#cccccc',
    borderwidth=1,
    borderpad=5,
    font=dict(size=18, family=FONT)
)

fig1.update_layout(
    title=dict(
        text="Actual vs Predicted: F<sub>ROH</sub> (%)",
        font=dict(size=18, family=FONT, color='#222222'),
        x=0.5, xanchor='center'
    ),
    xaxis=dict(
        title=dict(text="Actual F<sub>ROH</sub> (%)", font=dict(size=18, family=FONT)),
        rangemode='tozero',
        showgrid=True,
        gridcolor='#eeeeee',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    yaxis=dict(
        title=dict(text="Predicted F<sub>ROH</sub> (%)", font=dict(size=18, family=FONT)),
        rangemode='tozero',
        showgrid=True,
        gridcolor='#eeeeee',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    legend=dict(
        title=dict(text="<b>IUCN Status</b>", font=dict(size=18, family=FONT)),
        x=0.72, y=0.08,
        bgcolor='rgba(255,255,255,0.85)',
        bordercolor='#cccccc',
        borderwidth=1,
        font=dict(size=18, family=FONT)
    ),
    width=700, height=620,
    plot_bgcolor='white',
    paper_bgcolor='white',
    template='plotly_white',
    margin=dict(l=70, r=40, t=70, b=70),
    font=dict(family=FONT)
)


## -- Save figures -- ##
fig1.write_html(f"{output_dir}/HTML_Figures/{today}_actual_vs_predicted_FROH.html")
print(f"Saved: {output_dir}/HTML_Figures/{today}_actual_vs_predicted_FROH.html")

fig1.write_image(f"{output_dir}/Main_Figures/{today}_actual_vs_predicted_FROH.svg")
print(f"Saved: {output_dir}/Main_Figures/{today}_actual_vs_predicted_FROH.svg")
#### ---- End ---- ####


#### ---- Plot the coefficient strength for the FROH model ---- ####
params   = result.params.drop('const')
conf_int = result.conf_int().drop('const')
pvals        = result.pvalues.drop('const')

## -- Rename labesl for coefficients -- ##
display_names = [rename_coeff_label(n) for n in params.index]

n_params      = len(params)
sig_mask      = pvals < 0.05
point_colours = ['#1a6e8a' if s else '#b0b0b0' for s in sig_mask]
ci_colours    = ['#1a6e8a' if s else '#cccccc'  for s in sig_mask]

fig2 = go.Figure()

## -- Alternating background color -- ##
for i in range(n_params):
    fig2.add_shape(
        type='rect',
        xref='paper', yref='y',
        x0=0, x1=1,
        y0=i - 0.5, y1=i + 0.5,
        fillcolor='#f5f5f5' if i % 2 == 0 else 'white',
        line=dict(width=0),
        layer='below'
    )

## -- Add zero reference line -- ##
fig2.add_vline(
    x=0,
    line=dict(color='#444444', dash='dot', width=1.2)
)

## -- Add confidence intervals -- ##
for idx, (name, row) in enumerate(conf_int.iterrows()):
    col = ci_colours[idx]
    lw  = 3 if sig_mask.iloc[idx] else 1.8
    fig2.add_trace(go.Scatter(
        x=[row[0], row[1]],
        y=[idx, idx],
        mode='lines',
        line=dict(color=col, width=lw),
        showlegend=False,
        hoverinfo='skip'
    ))
    for xval in [row[0], row[1]]:
        fig2.add_trace(go.Scatter(
            x=[xval], y=[idx],
            mode='markers',
            marker=dict(color=col, size=4, symbol='line-ns-open',
                        line=dict(color=col, width=lw)),
            showlegend=False,
            hoverinfo='skip'
        ))

## -- Add coefficient points -- ##
fig2.add_trace(go.Scatter(
    x=params.values,
    y=list(range(n_params)),
    mode='markers',
    marker=dict(
        color=point_colours,
        size=10,
        line=dict(color='rgba(0,0,0,0.3)', width=0.8)
    ),
    text=[
        f"<b>{display_names[i]}</b><br>"
        f"Coef: {params.iloc[i]:.4f}<br>"
        f"95% CI: [{conf_int.iloc[i, 0]:.4f}, {conf_int.iloc[i, 1]:.4f}]<br>"
        f"p = {pvals.iloc[i]:.4f}"
        for i in range(n_params)
    ],
    hovertemplate="%{text}<extra></extra>",
    name='Coefficient',
    showlegend=False
))

## -- Add legend dummy traces -- ##
fig2.add_trace(go.Scatter(
    x=[None], y=[None], mode='markers',
    marker=dict(color='#1a6e8a', size=9),
    name='p < 0.05 (significant)'
))
fig2.add_trace(go.Scatter(
    x=[None], y=[None], mode='markers',
    marker=dict(color='#b0b0b0', size=9),
    name='p ≥ 0.05'
))

## -- Include R² annotation -- ##
fig2.add_annotation(
    xref='paper', yref='paper',
    x=0.98, y=0.01,
    text=f"<b>R² = {result.rsquared:.4f}</b>",
    showarrow=False,
    bgcolor='white',
    bordercolor='#cccccc',
    borderwidth=1,
    borderpad=5,
    font=dict(size=18, family=FONT),
    xanchor='right'
)

fig2.update_layout(
    title=dict(
        text="F<sub>ROH</sub> Model: Phylogenetically-Corrected Coefficients (PGLS)",
        font=dict(size=18, family=FONT, color='#222222'),
        x=0.5, xanchor='center'
    ),
    xaxis=dict(
        title=dict(
            text="Coefficient (proportional effect on F<sub>ROH</sub>)",
            font=dict(size=18, family=FONT)
        ),
        showgrid=True,
        gridcolor='#e0e0e0',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    yaxis=dict(
        tickmode='array',
        tickvals=list(range(n_params)),
        ticktext=display_names,
        tickfont=dict(size=11, family=FONT),
        showgrid=False
    ),
    height=max(650, n_params * 30),
    width=820,
    plot_bgcolor='white',
    paper_bgcolor='white',
    template='plotly_white',
    legend=dict(
        title=dict(text="<b>Significance</b>", font=dict(size=18, family=FONT)),
        x=0.75, y=0.99,
        bgcolor='rgba(255,255,255,0.85)',
        bordercolor='#cccccc',
        borderwidth=1,
        font=dict(size=11, family=FONT)
    ),
    margin=dict(l=200, r=40, t=70, b=60),
    font=dict(family=FONT)
)


# scaler_froh = StandardScaler()
# X_froh_scaled = pd.DataFrame(
#     scaler_froh.fit_transform(data_clean[fixed_effects].astype(float)),
#     columns=fixed_effects,
#     index=data_clean[fixed_effects].index
# )
# y_froh_scaled = pd.Series(
#     StandardScaler().fit_transform(y.values.reshape(-1, 1)).flatten(),
#     index=y.index
# )
# model_froh_scaled = sm.OLS(y_froh_scaled, sm.add_constant(X_froh_scaled)).fit()
# params_froh   = model_froh_scaled.params.drop('const')
# conf_int_froh = model_froh_scaled.conf_int().drop('const')
# pvals         = model_froh_scaled.pvalues.drop('const')

# ## -- Rename labels using the habitat map -- ##
# display_names = [rename_coeff_label(n) for n in params_froh.index]

# n_params      = len(params_froh)
# sig_mask      = pvals < 0.05
# point_colours = ['#1a6e8a' if s else '#b0b0b0' for s in sig_mask]
# ci_colours    = ['#1a6e8a' if s else '#cccccc'  for s in sig_mask]

# fig2 = go.Figure()

# ## -- Alternating background color -- ##
# for i in range(n_params):
#     fig2.add_shape(
#         type='rect',
#         xref='paper', yref='y',
#         x0=0, x1=1,
#         y0=i - 0.5, y1=i + 0.5,
#         fillcolor='#f5f5f5' if i % 2 == 0 else 'white',
#         line=dict(width=0),
#         layer='below'
#     )

# ## -- Add confidence intervals -- ##
# for idx, (name, row) in enumerate(conf_int_froh.iterrows()):
#     col = ci_colours[idx]
#     lw  = 3 if sig_mask.iloc[idx] else 1.8
#     fig2.add_trace(go.Scatter(
#         x=[row[0], row[1]],
#         y=[idx, idx],
#         mode='lines',
#         line=dict(color=col, width=lw),
#         showlegend=False,
#         hoverinfo='skip'
#     ))
#     for xval in [row[0], row[1]]:
#         fig2.add_trace(go.Scatter(
#             x=[xval], y=[idx],
#             mode='markers',
#             marker=dict(color=col, size=4, symbol='line-ns-open',
#                         line=dict(color=col, width=lw)),
#             showlegend=False,
#             hoverinfo='skip'
#         ))

# ## -- Add Coefficient points -- ##
# fig2.add_trace(go.Scatter(
#     x=params_froh.values,
#     y=list(range(n_params)),
#     mode='markers',
#     marker=dict(
#         color=point_colours,
#         size=10,
#         line=dict(color='rgba(0,0,0,0.3)', width=0.8)
#     ),
#     text=[
#         f"<b>{display_names[i]}</b><br>"
#         f"Coef: {params_froh.iloc[i]:.4f}<br>"
#         f"95% CI: [{conf_int_froh.iloc[i, 0]:.4f}, {conf_int_froh.iloc[i, 1]:.4f}]<br>"
#         f"p = {pvals.iloc[i]:.4f}"
#         for i in range(n_params)
#     ],
#     hovertemplate="%{text}<extra></extra>",
#     name='Coefficient',
#     showlegend=False
# ))

# ## -- Add zero reference line -- ##
# fig2.add_vline(
#     x=0,
#     line=dict(color='#444444', dash='dot', width=1.2)
# )

# ## -- Add legend dummy traces -- ##
# fig2.add_trace(go.Scatter(
#     x=[None], y=[None], mode='markers',
#     marker=dict(color='#1a6e8a', size=9),
#     name='p < 0.05 (significant)'
# ))
# fig2.add_trace(go.Scatter(
#     x=[None], y=[None], mode='markers',
#     marker=dict(color='#b0b0b0', size=9),
#     name='p ≥ 0.05'
# ))

# ## -- Include R² annotation -- ##
# fig2.add_annotation(
#     xref='paper', yref='paper',
#     x=0.98, y=0.01,
#     text=f"<b>R² = {model_froh_scaled.rsquared:.4f}</b>",
#     showarrow=False,
#     bgcolor='white',
#     bordercolor='#cccccc',
#     borderwidth=1,
#     borderpad=5,
#     font=dict(size=18, family=FONT),
#     xanchor='right'
# )

# fig2.update_layout(
#     title=dict(
#         text="F<sub>ROH</sub> Model: Standardised Coefficients",
#         font=dict(size=18, family=FONT, color='#222222'),
#         x=0.5, xanchor='center'
#     ),
#     xaxis=dict(
#         title=dict(
#             text="Standardised Coefficient (proportional effect on F<sub>ROH</sub>)",
#             font=dict(size=18, family=FONT)
#         ),
#         showgrid=True,
#         gridcolor='#e0e0e0',
#         zeroline=False,
#         tickfont=dict(size=11, family=FONT)
#     ),
#     yaxis=dict(
#         tickmode='array',
#         tickvals=list(range(n_params)),
#         ticktext=display_names,  
#         tickfont=dict(size=11, family=FONT),
#         showgrid=False
#     ),
#     height=max(650, n_params * 30),
#     width=820,
#     plot_bgcolor='white',
#     paper_bgcolor='white',
#     template='plotly_white',
#     legend=dict(
#         title=dict(text="<b>Significance</b>", font=dict(size=18, family=FONT)),
#         x=0.75, y=0.99,
#         bgcolor='rgba(255,255,255,0.85)',
#         bordercolor='#cccccc',
#         borderwidth=1,
#         font=dict(size=11, family=FONT)
#     ),
#     margin=dict(l=200, r=40, t=70, b=60),
#     font=dict(family=FONT)
# )

## -- Save figures -- ##
fig2.write_html(f"{output_dir}/HTML_Figures/{today}_regression_coeff_FROH.html")
print(f"Saved: {output_dir}/HTML_Figures/{today}_regression_coeff_FROH.html")

fig2.write_image(f"{output_dir}/Main_Figures/{today}_regression_coeff_FROH.svg")
print(f"Saved: {output_dir}/Main_Figures/{today}_regression_coeff_FROH.svg")
## -- END -- ##





########## ========== Get covariance and run Heterozygosity model ========== ##########
## -- Fit PGLS -- ##
fixed_effects = ['IUCN_encoded', 'all_aln_FROH_percent', 'Aqu', 'Fos', 'Ter', 'Arb', 'Aer'] + hab_cols
X = sm.add_constant(data_clean[fixed_effects].astype(float))
y = data_clean['heterozygosity'].astype(float)
y_logit = pd.Series(
    logit(y.values / 100).flatten(),
    index=y.index
)

het_model = GLS(y_logit, X, sigma=cov_array)
het_result = het_model.fit()
print(result.summary())


########## ========== Plot heterozygosity model results ========== ##########
#### ---- Plot actual vs. predicted for heterozygosity model ---- ####
y_het_pred   = expit(result.fittedvalues) * 100
y_het_actual = y.values

min_val = min(y_het_actual.min(), y_het_pred.min())
max_val = max(y_het_actual.max(), y_het_pred.max())

fig1 = go.Figure()

## -- Plot 1:1 reference line -- ##
fig1.add_trace(go.Scatter(
    x=[min_val, max_val],
    y=[min_val, max_val],
    mode='lines',
    line=dict(color='#888888', dash='dash', width=1.5),
    name='1 : 1 line',
    hoverinfo='skip'
))

## -- Add IUCN traces -- ##
for status in iucn_order:
    mask = data_clean['IUCN'] == status
    if not mask.any():
        continue

    subset     = data_clean[mask]
    actual_sub = y_het_actual[mask]
    pred_sub   = y_het_pred[mask]

    hover_sub = [
        f"<b>{row['species']}</b><br>"
        f"IUCN: {row['IUCN']}<br>"
        f"Lineage: {row['Extended_lineage']}<br>"
        f"FROH: {row['all_aln_FROH_percent']:.2f}%"
        f"Heterozygosity: {row['heterozygosity']}"
        for _, row in subset.iterrows()
    ]

    fig1.add_trace(go.Scatter(
        x=actual_sub,
        y=pred_sub,
        mode='markers',
        marker=dict(
            color=iucn_colors[status],
            symbol=iucn_symbols[status],
            size=10,
            opacity=0.9,
            line=dict(color='rgba(0,0,0,0.45)', width=0.8)
        ),
        text=hover_sub,
        hovertemplate="%{text}<br><b>Actual:</b> %{x:.2f}%<br><b>Predicted:</b> %{y:.2f}%<extra></extra>",
        name=status
    ))

## -- R² annotation -- ##
fig1.add_annotation(
    xref='paper', yref='paper',
    x=0.05, y=0.95,
    text=f"<b>R² = {het_result.rsquared:.4f}</b>",
    showarrow=False,
    bgcolor='white',
    bordercolor='#cccccc',
    borderwidth=1,
    borderpad=5,
    font=dict(size=18, family=FONT)
)

fig1.update_layout(
    title=dict(
        text="Actual vs Predicted: Heterozygosity per kb outside of ROH",
        font=dict(size=18, family=FONT, color='#222222'),
        x=0.5, xanchor='center'
    ),
    xaxis=dict(
        title=dict(text="Actual Heterozygosity per kb outside of ROH", font=dict(size=18, family=FONT)),
        range=[0, 10], 
        showgrid=True,
        gridcolor='#eeeeee',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    yaxis=dict(
        title=dict(text="Predicted heterozygosity per kb outside of ROH", font=dict(size=18, family=FONT)),
        showgrid=True,
        gridcolor='#eeeeee',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    legend=dict(
        title=dict(text="<b>IUCN Status</b>", font=dict(size=18, family=FONT)),
        x=0.72, y=0.08,
        bgcolor='rgba(255,255,255,0.85)',
        bordercolor='#cccccc',
        borderwidth=1,
        font=dict(size=18, family=FONT)
    ),
    width=700, height=620,
    plot_bgcolor='white',
    paper_bgcolor='white',
    template='plotly_white',
    margin=dict(l=70, r=40, t=70, b=70),
    font=dict(family=FONT)
)


## -- Save figures -- ##
fig1.write_html(f"{output_dir}/HTML_Figures/{today}_actual_vs_predicted_Het.html")
print(f"Saved: {output_dir}/HTML_Figures/{today}_actual_vs_predicted_Het.html")

fig1.write_image(f"{output_dir}/Supplementary_Figures/{today}_actual_vs_predicted_Het.svg")
print(f"Saved: {output_dir}/Supplementary_Figures/{today}_actual_vs_predicted_Het.svg")
#### ---- End ---- ####


#### ---- Plot the coefficient strength for the Heterozygosity model ---- ####
params_het   = het_result.params.drop('const')
conf_int_het = het_result.conf_int().drop('const')
pvals        = het_result.pvalues.drop('const')

## -- Rename labesl for coefficients -- ##
display_names = [rename_coeff_label(n) for n in params_het.index]

n_params      = len(params_het)
sig_mask      = pvals < 0.05
point_colours = ['#1a6e8a' if s else '#b0b0b0' for s in sig_mask]
ci_colours    = ['#1a6e8a' if s else '#cccccc'  for s in sig_mask]

fig2 = go.Figure()

## -- Alternating background color -- ##
for i in range(n_params):
    fig2.add_shape(
        type='rect',
        xref='paper', yref='y',
        x0=0, x1=1,
        y0=i - 0.5, y1=i + 0.5,
        fillcolor='#f5f5f5' if i % 2 == 0 else 'white',
        line=dict(width=0),
        layer='below'
    )

## -- Add zero reference line -- ##
fig2.add_vline(
    x=0,
    line=dict(color='#444444', dash='dot', width=1.2)
)

## -- Add confidence intervals -- ##
for idx, (name, row) in enumerate(conf_int_het.iterrows()):
    col = ci_colours[idx]
    lw  = 3 if sig_mask.iloc[idx] else 1.8
    fig2.add_trace(go.Scatter(
        x=[row[0], row[1]],
        y=[idx, idx],
        mode='lines',
        line=dict(color=col, width=lw),
        showlegend=False,
        hoverinfo='skip'
    ))
    for xval in [row[0], row[1]]:
        fig2.add_trace(go.Scatter(
            x=[xval], y=[idx],
            mode='markers',
            marker=dict(color=col, size=4, symbol='line-ns-open',
                        line=dict(color=col, width=lw)),
            showlegend=False,
            hoverinfo='skip'
        ))

## -- Add coefficient points -- ##
fig2.add_trace(go.Scatter(
    x=params_het.values,
    y=list(range(n_params)),
    mode='markers',
    marker=dict(
        color=point_colours,
        size=10,
        line=dict(color='rgba(0,0,0,0.3)', width=0.8)
    ),
    text=[
        f"<b>{display_names[i]}</b><br>"
        f"Coef: {params_het.iloc[i]:.4f}<br>"
        f"95% CI: [{conf_int_het.iloc[i, 0]:.4f}, {conf_int_het.iloc[i, 1]:.4f}]<br>"
        f"p = {pvals.iloc[i]:.4f}"
        for i in range(n_params)
    ],
    hovertemplate="%{text}<extra></extra>",
    name='Coefficient',
    showlegend=False
))

## -- Add legend dummy traces -- ##
fig2.add_trace(go.Scatter(
    x=[None], y=[None], mode='markers',
    marker=dict(color='#1a6e8a', size=9),
    name='p < 0.05 (significant)'
))
fig2.add_trace(go.Scatter(
    x=[None], y=[None], mode='markers',
    marker=dict(color='#b0b0b0', size=9),
    name='p ≥ 0.05'
))

## -- Include R² annotation -- ##
fig2.add_annotation(
    xref='paper', yref='paper',
    x=0.98, y=0.01,
    text=f"<b>R² = {het_result.rsquared:.4f}</b>",
    showarrow=False,
    bgcolor='white',
    bordercolor='#cccccc',
    borderwidth=1,
    borderpad=5,
    font=dict(size=18, family=FONT),
    xanchor='right'
)

fig2.update_layout(
    title=dict(
        text="Heterozygosity Model: Phylogenetically-Corrected Coefficients (PGLS)",
        font=dict(size=18, family=FONT, color='#222222'),
        x=0.5, xanchor='center'
    ),
    xaxis=dict(
        title=dict(
            text="Coefficient (effect on logit heterozygosity outside of ROH)",
            font=dict(size=18, family=FONT)
        ),
        showgrid=True,
        gridcolor='#e0e0e0',
        zeroline=False,
        tickfont=dict(size=11, family=FONT)
    ),
    yaxis=dict(
        tickmode='array',
        tickvals=list(range(n_params)),
        ticktext=display_names,
        tickfont=dict(size=11, family=FONT),
        showgrid=False
    ),
    height=max(650, n_params * 30),
    width=820,
    plot_bgcolor='white',
    paper_bgcolor='white',
    template='plotly_white',
    legend=dict(
        title=dict(text="<b>Significance</b>", font=dict(size=18, family=FONT)),
        x=0.75, y=0.99,
        bgcolor='rgba(255,255,255,0.85)',
        bordercolor='#cccccc',
        borderwidth=1,
        font=dict(size=11, family=FONT)
    ),
    margin=dict(l=200, r=40, t=70, b=60),
    font=dict(family=FONT)
)


## -- Save figures -- ##
fig2.write_html(f"{output_dir}/HTML_Figures/{today}_regression_coeff_Het.html")
print(f"Saved: {output_dir}/HTML_Figures/{today}_regression_coeff_Het.html")

fig2.write_image(f"{output_dir}/Supplementary_Figures/{today}_regression_coeff_Het.svg")
print(f"Saved: {output_dir}/Supplementary_Figures/{today}_regression_coeff_Het.svg")
## -- END -- ##
