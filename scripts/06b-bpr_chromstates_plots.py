#!/usr/bin/env python

## PURPOSE: Plot chrstates.
## USAGE:   python 06b-bpr_chromstates_plots.py <chrstates>
## OUTPUT:  chromstates.pdf

## LIBRARIES
import os
os.environ['MPLCONFIGDIR']= "/projects/lansdorp/lansdorp_scratch/nfwork/$USER"

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Plots
chrstates = pd.read_table(sys.argv[1]).replace('?', "Unknown").fillna('Unknown')
chrOrder = chrstates.columns.tolist()[1:]
chrstates["Sample"] = chrstates['library'].str.split('-r\d+-c\d+').str[0]

melt_chrstates = chrstates.drop(columns=['library']).melt(id_vars='Sample', var_name='Chromosome', value_name='State')
count_df = pd.pivot_table(melt_chrstates, values='State', index=['Sample','Chromosome'], columns='State', aggfunc='size')[['wc', 'ww', 'cc', 'Unknown']]
row_sums = count_df.sum(axis=1)
pct_df = (count_df.div(row_sums, axis=0) * 100).reset_index()

pdf = PdfPages('chromstates.pdf')

with PdfPages('chromstates.pdf') as pdf:
    for sample in set(pct_df["Sample"]) - set(["empty", "negativecontrol"]):
        sub_df = pct_df.query('Sample == @sample').set_index("Chromosome").reindex(index=chrOrder).reset_index()
        
        sub_df.plot(kind='bar', x='Chromosome', stacked=True, figsize=(15,10), color=['firebrick', 'steelblue', 'darkorange', 'grey'])
        plt.xlabel('Chromosome', fontsize=20)
        plt.ylabel('Percentage of States', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.legend(bbox_to_anchor=(1.01,1), loc='upper left', fontsize=15)
        plt.title(sample, fontsize=25)
        plt.tight_layout()
        pdf.savefig()
        plt.close()