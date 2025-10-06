#!/usr/bin/env python

## PURPOSE: Plot heatmaps from metric_details.tsv
## USAGE:   python 10b-metrics_plots.py <metrics_file>
## OUTPUT:  heatmaps.pdf

## LIBRARIES
import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

## VARIABLE
metrics_file = sys.argv[1]

## SCRIPT
dim = 15

columns_order = ['Quality', 'Background', 'Duplication_rate', 'Complexity_at_1Gb', 'Median_insert_size', 
                 'Adapters_at_least_50bp', 'Alignment_rate', 'Initial_reads_aligned', 'Processed_reads_aligned', 
                 'Coverage', 'Reads_per_Mb', 'Mean_GC', 'Percent_WC', 'Read_length']

# Data
metdf = pd.read_table(metrics_file)
metdf_chiponly = metdf.query("not R.isna()").copy()
metdf_chiponly[['R', 'C']] = metdf_chiponly[['R', 'C']].astype({'R': int, 'C': int})
metdf_chiponly['RC'] = metdf_chiponly.apply(lambda row: [row['R'], row['C']], axis=1)
metdf_sampleonly = metdf_chiponly[~metdf_chiponly['Sample'].str.contains('negative|empty', case=False, na=False)][['Sample', 'Quality']]
# metdf_sampleonly = metdf_chiponly.query('Sample != ["negativecontrol", "empty"]')[['Sample', 'Quality']]

# Plots
pdf = PdfPages('heatmaps.pdf')

# Sample Map
# labels = metdf_chiponly.groupby('Sample').agg(R=('R', 'min'), C=('C','min'), N=('Sample', 'count')).reset_index().reset_index()
labels = metdf_chiponly.groupby('Sample').agg(RC=('RC', 'min'), N=('Sample', 'count')).reset_index().reset_index()
labels[['R', 'C']] = pd.DataFrame(labels['RC'].to_list())
labels['Label'] = labels['Sample']  + " (" + labels['N'].astype(str) +")"
try:
    labels.loc[labels['Sample'].str.contains("negative|empty"), 'index'] = -1
except IndexError:
    pass

samplemap = metdf_chiponly[['Sample', 'R', 'C']].merge(labels[['Sample', 'index']]).pivot(index="R", columns="C", values="index")
combined_cmap = colors.ListedColormap(np.vstack((plt.cm.tab20(np.linspace(0, 1, 20)), plt.cm.tab20b(np.linspace(0, 1, 20)), plt.cm.Set1(np.linspace(0, 1, 9)), plt.cm.Dark2(np.linspace(0, 1, 8)))))

allcolumns = list(range(metdf_chiponly['C'].min(), metdf_chiponly['C'].max()+1))
allrows = list(range(metdf_chiponly['R'].min(), metdf_chiponly['R'].max()+1))
samplemap = samplemap.reindex(columns=allcolumns).reindex(allrows)

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
        matrix = matrix.reindex(columns=allcolumns).reindex(allrows)
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
