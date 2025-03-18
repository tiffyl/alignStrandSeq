#!/usr/bin/env python

## PURPOSE: Plot background maps.
## USAGE:   python 09-background-plot.py <paired>
## OUTPUT:  background.pdf

## LIBRARIES
import os
import sys
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

## VARIABLES
paired = bool(sys.argv[1])

## SCRIPT
def density_plot(name):
    background = pd.read_table(name + "background.stats")
    directional = pd.read_table(name + "directional.stats")
    readcounts = pd.read_table(name + "readcounts.txt", sep=" ").columns
    
    both=pd.concat([background, directional])
    
    fig, ax = plt.subplots(1,2)
    
    # Heatmap plots for the two panels
    sns.kdeplot(data=background, x=background.columns[0], y=background.columns[1], cmap="Reds", fill=True, bw_adjust=.5, ax=ax[0]).set(xlabel=None, ylabel=None)
    sns.kdeplot(data=directional, x=directional.columns[0], y=directional.columns[1], cmap="Reds", fill=True, bw_adjust=.5, ax=ax[1]).set(xlabel=None, ylabel=None)
    
    # Want panels to have same scale
    plt.setp(ax, xlim=(0, 1), ylim=(both.min()["12_seq_len"], both.max()["12_seq_len"]))
    
    ax[0].title.set_text('Background')
    ax[1].title.set_text('Directional')
    ax[1].get_yaxis().set_visible(False)
    
    if name != "":
        name = "\n" + name.strip(".")
    
    fig.suptitle("Comparison of background and directional reads on chr1q%s" % name)
    fig.supxlabel("GC content")
    
    # Paired end reads---we can't look at insert size there
    if paired:
        fig.supylabel("Insert size (bp)")
    else:
        fig.supylabel("Read length (bp)")
    
    # Annotating the plots with fragment counts
    percent = round(int(readcounts[0]) / int(readcounts[1]) * 100, 1)
    text_background= readcounts[0] + " ("+ str(percent) + "%)"
    text_directional= readcounts[1] + ' fragments'
    ax[0].text(0.05,0.03, text_background, transform=ax[0].transAxes)
    ax[1].text(0.05,0.03, text_directional, transform=ax[1].transAxes)

    plt.tight_layout()
    pdf.savefig()
    plt.close()

samples = sorted(set([ os.path.basename(x).split(".")[0] for x in glob.glob("*.*.stats") ]))

with PdfPages('background.pdf') as pdf:
    if os.path.isfile("background.stats") & os.path.isfile("directional.stats"):
        density_plot("")
        
    for s in samples:
        if os.path.isfile(s + ".background.stats") & os.path.isfile(s + ".directional.stats"):
            density_plot(s + ".")
