#### ==== Script to calculate Linear Regression Model to predict ROH ==== ####
## Date: 20260527 (May 27th, 2026)
## Author: Amanda Gardiner
## Version 3 (V2 is 20260327_ROH_Linear_Regression_V2.py)
## GOAL: Create a linear regression model to see what predicts ROH for a species
## NOTES:
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
import re as re
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
    """Extract first number from each comma-separated habitat code"""
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


########## ========== Get covariance and run model ========== ##########
## -- Match species in data to tips in tree -- ##
def build_phylo_cov_matrix(tree, accession_list):
    """
    Build a covariance matrix based on shared branch lengths.
    accession_list must match tree tip names exactly.
    """
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


########## ========== Plot results ========== ##########
#### ---- Plot actual vs. predicted for FROH model ---- ####
y_froh_pred = expit(result.fittedvalues) * 100
y_froh_actual = y.values

hover_text = [
    f"<b>{row['species']}</b><br>"
    f"IUCN: {row['IUCN']}<br>"
    f"Lineage: {row['Extended_lineage']}<br>"
    f"FROH: {row['all_aln_FROH_percent']:.2f}%"
    for _, row in data_clean.iterrows()
]

min_val = min(y_froh_actual.min(), y_froh_pred.min())
max_val = max(y_froh_actual.max(), y_froh_pred.max())

fig1 = go.Figure()
fig1.add_trace(go.Scatter(
    x=y_froh_actual,
    y=y_froh_pred,
    mode='markers',
    marker=dict(
        color='lightblue',
        line=dict(color='steelblue', width=0.8),
        size=8,
        opacity=0.8
    ),
    text=hover_text,
    hovertemplate="%{text}<br><b>Actual:</b> %{x:.2f}%<br><b>Predicted:</b> %{y:.2f}%<extra></extra>",
    name='Samples'
))
fig1.add_trace(go.Scatter(
    x=[min_val, max_val],
    y=[min_val, max_val],
    mode='lines',
    line=dict(color='red', dash='dash', width=1.5),
    name='Perfect fit',
    hoverinfo='skip'
))
fig1.add_annotation(
    xref='paper', yref='paper',
    x=0.05, y=0.95,
    text=f"R² = {result.rsquared:.4f}",
    showarrow=False,
    bgcolor='white',
    bordercolor='grey',
    borderwidth=1,
    font=dict(size=12)
)
fig1.update_layout(
    title="Actual vs Predicted: all_aln_FROH_percent",
    xaxis_title="Actual FROH (%)",
    yaxis_title="Predicted FROH (%)",
    xaxis=dict(rangemode='tozero'),
    yaxis=dict(rangemode='tozero'),
    width=700, height=600,
    legend=dict(x=0.7, y=0.05),
    template='plotly_white'
)
fig1.write_html(f"Figures/{today}_actual_vs_predicted_FROH.html")
print(f"Saved: Figures/{today}_actual_vs_predicted_FROH.html")
#### ---- End ---- ####


#### ---- Plot the coefficient strength for the FROH model ---- ####
scaler_froh = StandardScaler()
X_froh_scaled = pd.DataFrame(
    scaler_froh.fit_transform(data_clean[fixed_effects].astype(float)),
    columns=fixed_effects,
    index=data_clean[fixed_effects].index
)
y_froh_scaled = pd.Series(
    StandardScaler().fit_transform(y.values.reshape(-1, 1)).flatten(),
    index=y.index
)
model_froh_scaled = sm.OLS(y_froh_scaled, sm.add_constant(X_froh_scaled)).fit()
params_froh = model_froh_scaled.params.drop('const')
conf_int_froh = model_froh_scaled.conf_int().drop('const')

## -- Colour points by significance -- ##
pvals = model_froh_scaled.pvalues.drop('const')
point_colours = ['#1a9988' if p < 0.05 else '#aaaaaa' for p in pvals]

fig2 = go.Figure()

## -- Confidence interval lines -- ##
for idx, (name, row) in enumerate(conf_int_froh.iterrows()):
    fig2.add_trace(go.Scatter(
        x=[row[0], row[1]],
        y=[idx, idx],
        mode='lines',
        line=dict(color='#eb5600', width=3),
        showlegend=False,
        hoverinfo='skip'
    ))

## -- Coefficient points -- ##
fig2.add_trace(go.Scatter(
    x=params_froh.values,
    y=list(range(len(params_froh))),
    mode='markers',
    marker=dict(color=point_colours, size=8),
    text=[
        f"<b>{name}</b><br>"
        f"Coef: {params_froh[name]:.4f}<br>"
        f"95% CI: [{conf_int_froh.loc[name, 0]:.4f}, {conf_int_froh.loc[name, 1]:.4f}]<br>"
        f"p = {pvals[name]:.4f}"
        for name in params_froh.index
    ],
    hovertemplate="%{text}<extra></extra>",
    name='Coefficient',
    showlegend=False
))

## -- Zero line -- ##
fig2.add_vline(x=0, line=dict(color='#eb5600', dash='dash'))

## -- Dummy traces for legend -- ##
fig2.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
    marker=dict(color='#1a9988', size=8), name='p < 0.05'))
fig2.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
    marker=dict(color='#aaaaaa', size=8), name='p ≥ 0.05'))

fig2.add_annotation(
    xref='paper', yref='paper',
    x=0.98, y=0.02,
    text=f"R² = {model_froh_scaled.rsquared:.4f}",
    showarrow=False,
    bgcolor='white',
    bordercolor='grey',
    borderwidth=1,
    font=dict(size=11),
    xanchor='right'
)
fig2.update_layout(
    title="FROH Model: Strength of Relationships",
    xaxis_title="Standardised Coefficient (Proportional Effect)",
    yaxis=dict(
        tickmode='array',
        tickvals=list(range(len(params_froh))),
        ticktext=list(params_froh.index),
        tickfont=dict(size=9)
    ),
    height=max(600, len(params_froh) * 28),
    width=750,
    template='plotly_white',
    legend=dict(x=0.75, y=0.99)
)
fig2.write_html(f"Figures/{today}_regression_coeff_FROH.html")
print(f"Saved: Figures/{today}_regression_coeff_FROH.html")
## -- END -- ##










