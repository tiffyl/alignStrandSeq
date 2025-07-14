#!/usr/bin/env python

## PURPOSE: Plot preseq estimates as a curve for libraries with preseq.estimate in current directory. 
## USAGE:   python 07b-preseq-sample-plot.py <genomesize> <preseqestdir> 
## OUTPUT:  {sampleId}.preseq.pdf preseq_all.pdf 

## LIBRARIES
import sys
import os
import re
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 

## VARIABLES
genomesize = float(sys.argv[1])
preseqdir = sys.argv[2]

limit = 2
step = 0.2
dim=12

## SCRIPT
# Read all preseq estimates into a single file
preseqlist = glob.glob(preseqdir + "/*.preseq.estimate")

allpreseqdf = pd.DataFrame()

for file in preseqlist:
    libId = os.path.basename(file).split(".")[0]
    df = pd.read_table(file)
    df["Library"] = libId

    # logfile = os.path.dirname(file) + "/preseq." + libId + ".log"
    logfile = preseqdir + libId + ".preseq.log"
    with open(logfile) as log:
        for line in log:
            if re.search("TOTAL COVERED BASES", line):
                df["Breadth"] = float(line.strip().split("=")[1].strip())
                
    allpreseqdf = pd.concat([allpreseqdf, df]).reset_index(drop=True)

allpreseqdf['Total_Gb'] = allpreseqdf['TOTAL_BASES'] / 1000000000
allpreseqdf['Genome_Coverage'] = allpreseqdf['EXPECTED_COVERED_BASES'] / genomesize * 100
allpreseqdf['Lower'] = allpreseqdf['LOWER_95%CI'] / genomesize * 100
allpreseqdf['Upper'] = allpreseqdf['UPPER_95%CI'] / genomesize * 100
allpreseqdf["SampleId"] = allpreseqdf['Library'].str.split('-r..-c..', expand=True)[0]
allpreseqdf = allpreseqdf.iloc[:, 4:]

# Samples
samples = allpreseqdf["SampleId"].unique().astype(str)
samples = samples[~np.char.find(samples, 'neg') >= 0]

# Plot
def complexity_plot(sampleId):
    if (sampleId == ""):
        preseqdf = allpreseqdf
        outfile = "preseq_all.pdf"
        
        fig = plt.figure(figsize=(dim + (dim/3), dim))
        gs = GridSpec(1,2, width_ratios=[3,1])
        ax = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
    else:
        preseqdf = allpreseqdf.query("SampleId == @sampleId")
        outfile = sampleId + ".preseq.pdf"
        
        fig, ax = plt.subplots(figsize=(dim, dim))

    max_gencov = max(preseqdf['Genome_Coverage'])
    breadth = preseqdf[["Library", "Breadth"]].drop_duplicates()["Breadth"]
    mean_breadth = np.mean(breadth) / 1000000000
    std_breadth = np.std(breadth) / 1000000000
    
    mean_preseqdf = preseqdf.iloc[:, 2:4].groupby("Total_Gb").mean().reset_index()
    mean_preseqdf['stdev'] = preseqdf.iloc[:, 2:4].groupby("Total_Gb").std()["Genome_Coverage"].tolist()
    mean_preseqdf['Lower'] = mean_preseqdf['Genome_Coverage'] - mean_preseqdf['stdev']
    mean_preseqdf['Upper'] = mean_preseqdf['Genome_Coverage'] + mean_preseqdf['stdev']
        
    for lib, libdf in preseqdf.groupby('Library'):
        ax.plot(libdf['Total_Gb'], libdf['Genome_Coverage'], color="lightgrey", label=lib, zorder=1)

    ax.plot(mean_preseqdf['Total_Gb'], mean_preseqdf['Genome_Coverage'], color="cornflowerblue", zorder=2)
    ax.fill_between(mean_preseqdf["Total_Gb"], mean_preseqdf["Upper"], mean_preseqdf["Lower"], color="cornflowerblue", alpha=0.15, zorder=2)    

    ax.plot(pd.Series([mean_breadth, mean_breadth]), pd.Series([0, max_gencov + 1]), linestyle='--', color="black")
    ax.fill_between(pd.Series([mean_breadth - std_breadth, mean_breadth + std_breadth]), 
                    pd.Series([0, 0]), 
                    pd.Series([max_gencov + 1, max_gencov + 1]), 
                    color="orange", alpha=0.25)
    ax.text(mean_breadth + std_breadth + 0.01, max_gencov - 1, f"Mean Observed Sequencing Effort: {round(mean_breadth, 2)}", fontsize=12)

    # Axis
    ax.set_xlim(0, limit)
    ax.set_xticks(np.arange(0, limit + step, step))
    ax.set_ylim(0, max_gencov)
    ax.set_yticks(np.arange(0, max_gencov + 1, 2))
    
    # Labels
    ax.set_title(sampleId, fontsize=20)
    ax.set_xlabel("Total Sequencing Effort (Gb)", fontsize=15)
    ax.set_ylabel("Predicted Breadth of Coverage (%)", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Histogram 
    if (sampleId == ""):
        hist1gb = preseqdf.query('Total_Gb == 1.0')[['Genome_Coverage']]
        
        ax2.hist(hist1gb, bins=list(range(int(round(max_gencov, 0)))), orientation="horizontal")
        ax2.set_ylim(0, max_gencov)
        ax2.set_yticks(np.arange(0, max_gencov + 1, 2))
            
        ax2.set_title("1GB Coverage", fontsize=20)
        ax2.set_xlabel("Number of Libraries", fontsize=15)
        ax2.tick_params(axis='both', labelsize=12)

    plt.tight_layout()
    plt.savefig(outfile, format='pdf', bbox_inches='tight')

# Per Sample
for sample in samples:
    complexity_plot(sample)

# All Libraries
complexity_plot("")