#!/usr/bin/env python

## PURPOSE: Group individual metrics_* files into a spreadsheet, and plot heatmaps
## USAGE:   python 10a-metrics_plots.py
## OUTPUT:  heatmaps.pdf

## LIBRARIES
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

## SCRIPT
dim = 15

columns_order = ['Quality', 'Background', 'Duplication_rate', 'Complexity_at_1Gb', 'Median_insert_size', 
                 'Adapters_at_least_50bp', 'Alignment_rate', 'Initial_reads_aligned', 'Processed_reads_aligned', 
                 'Coverage', 'Reads_per_Mb', 'Mean_GC', 'Percent_WC', 'Read_length']

# Data
metdf = pd.read_table('metrics_details.tsv')
metdf_chiponly = metdf.query("not R.isna()").copy()
metdf_chiponly[['R', 'C']] = metdf_chiponly[['R', 'C']].astype({'R': int, 'C': int})
metdf_sampleonly = metdf_chiponly.query('Sample != ["negativecontrol", "empty"]')[['Sample', 'Quality']]

# Plots
pdf = PdfPages('heatmaps.pdf')

# Sample Map
labels = metdf_chiponly.groupby('Sample').agg(R=('R', 'min'), C=('C','min'), N=('Sample', 'count')).reset_index().reset_index()
labels['Label'] = labels['Sample']  + " (" + labels['N'].astype(str) +")"
try:
    labels.loc[labels.query('Sample == "negativecontrol"').index[0], 'index'] = -1
    labels.loc[labels.query('Sample == "empty"').index[0], 'index'] = -1
except IndexError:
    pass

samplemap = metdf_chiponly[['Sample', 'R', 'C']].merge(labels[['Sample', 'index']]).pivot(index="R", columns="C", values="index")
combined_cmap = colors.ListedColormap(np.vstack((plt.cm.tab20(np.linspace(0, 1, 20)), plt.cm.tab20b(np.linspace(0, 1, 20)), plt.cm.Set1(np.linspace(0, 1, 9)), plt.cm.Dark2(np.linspace(0, 1, 8)))))

# Metrics Heatmaps
with PdfPages('heatmaps.pdf') as pdf:
    plt.figure(figsize=(dim, dim))
    sns.heatmap(samplemap, mask=(samplemap < 0), cmap=colors.ListedColormap(combined_cmap(np.linspace(0, 1, len(labels)))), cbar=False)
    sns.heatmap(samplemap, mask=(samplemap > -1), cmap='binary', vmin=-2, vmax=-1, cbar=False)

    for i in labels.query('index > -1').index.tolist():
        plt.text(labels.iloc[i]['C'] - min(labels['C']) + 0.5, 
                 labels.iloc[i]['R'] - min(labels['R']) + 0.5, 
                 labels.iloc[i]['Label'], fontsize=15, ha='left')

    plt.xlabel("Column", fontsize=20)
    plt.ylabel("Row", fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title("Sample Map", fontsize=25)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    for col in columns_order:
        plt.figure(figsize=(dim + 3, dim))
        matrix = metdf_chiponly.merge(labels[['Sample', 'index']]).pivot(index="R", columns="C", values=col)
        ax = sns.heatmap(matrix, mask=(matrix.isna()), cmap="viridis", vmin=0)
        
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=15)
        
        plt.xlabel("Column", fontsize=20)
        plt.ylabel("Row", fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.title(col, fontsize=25)

        plt.tight_layout()
        pdf.savefig()
        plt.close()
