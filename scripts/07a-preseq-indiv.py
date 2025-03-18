#!/usr/bin/env python

## PURPOSE: Plot preseq estimates as a curve per single cell.
## USAGE:   python 07a-preseq-plot.py <estimatefile> <logfile> <genomesize>
## OUTPUT:  {libId}.preseq.pdf {libId}.complexity.txt

## LIBRARIES
import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## VARIABLES
estfile = sys.argv[1]
logfile = sys.argv[2]
genomesize = float(sys.argv[3])

libId = os.path.basename(estfile).split(".")[0]

limit = 2
step = 0.2
dim=12

## SCRIPT
# Extract depth of library
with open(logfile) as file:
    for line in file:
        if re.search("TOTAL COVERED BASES", line):
            coveredbases = float(line.strip().split("=")[1].strip())

# Read table to convert to % for plotting
preseq = pd.read_table(estfile)
preseq['Total_Gb'] = preseq['TOTAL_BASES'] / 1000000000
preseq['Genome_Coverage'] = preseq['EXPECTED_COVERED_BASES'] / genomesize * 100
preseq['Lower'] = preseq['LOWER_95%CI'] / genomesize * 100
preseq['Upper'] = preseq['UPPER_95%CI'] / genomesize * 100

genfrac = coveredbases / genomesize
point_estimate = preseq.query('Total_Gb == 1.0')['Genome_Coverage'].values[0]/100

with open(libId + ".complexity.txt", "w") as outfile:
    outfile.write(f"{libId}\t{point_estimate}\n")

# Benchmark for 200 IGSR libraries (Chaisson et al. 2019)
benchmarkfile = pd.read_table('/projects/lansdorp/nextflow_pipelines/alignStrandSeq/scripts/files/preseq_benchmark.plot').dropna(axis=1)
std_benchmark = benchmarkfile.iloc[:, 1:].std(axis=1)
mean_benchmark = benchmarkfile.iloc[:, 1:].mean(axis=1)
benchmark = pd.DataFrame({"Total_Gb": benchmarkfile.iloc[:, 0], 
                          "Genome_Coverage": mean_benchmark, 
                          "Lower": mean_benchmark - std_benchmark, 
                          "Upper": mean_benchmark + std_benchmark})

# Breadth
breadthdf = pd.DataFrame({'Total_Gb':preseq['Total_Gb'], 'Breadth':genfrac * 100})
breadth5per = round(np.interp(5, preseq['Genome_Coverage'], preseq['Total_Gb']), 3)

# Plot
def matplot_complexity_curves(ax, df, col):    
    ax.plot(df["Total_Gb"], df["Genome_Coverage"], color=col)
    ax.fill_between(df["Total_Gb"], df["Upper"], df["Lower"], color=col, alpha=0.25)

fig, ax = plt.subplots(figsize=(dim, dim))
matplot_complexity_curves(ax, preseq, "cornflowerblue")
matplot_complexity_curves(ax, benchmark, "grey")

# Observed Breadth of Coverage
ax.plot(breadthdf['Total_Gb'], breadthdf['Breadth'], label="Observed Breadth", color='black', linestyle='--', alpha=0.5)
ax.text(breadthdf.iloc[10,0], breadthdf.iloc[1,1], f"Observed Breadth of Coverage: {round(breadthdf.iloc[1,1], 2)}%", 
        fontsize=15, verticalalignment='bottom', horizontalalignment='left', color='black')

# Axis
max_gencov = max(preseq['Genome_Coverage'])

ax.set_xlim(0, limit)
ax.set_xticks(np.arange(0, limit + step, step))
ax.set_ylim(0, max_gencov + 1)
ax.set_yticks(np.arange(0, max_gencov + 1, 2))
        
# Labels
ax.set_title(libId, fontsize=20)
ax.set_xlabel("Total Sequencing Effort (Gb)", fontsize=15)
ax.set_ylabel("Predicted Breadth of Coverage (%)", fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=12)

# Caption
caption = (f"Complexity curve and 95% CI (blue). This library requires {breadth5per}Gb sequencing for 5% breadth of coverage.\n"
           "The grey represents the mean and standard deviation of the HGSVC benchmark dataset \n"
           "(200 arbitrarily chosen libraries from Chaisson et al. 2019)")
plt.figtext(0.12, -0.01, caption, ha="left", fontsize=12, wrap=True)

plt.savefig(libId + ".preseq.pdf", format='pdf', bbox_inches='tight')

#####################
# Altair Plot
# import altair as alt
# dimension=800

# def complexity_curves(df, col):    
#     df = df.query("Total_Gb <= @limit")
    
#     chart = (
#         alt.Chart(df).mark_line(color=col).encode(
#             alt.X("Total_Gb", scale=alt.Scale(domainMax=limit, clamp=True), 
#                   axis=alt.Axis(values=np.arange(0, limit + step, step)),
#                   title="Total Sequencing Effort (Gb)"),
#             alt.Y("Genome_Coverage", 
#                   title="Predicted Breadth of Coverage (%)", 
#                   axis=alt.Axis(values=np.arange(0, max(df['Genome_Coverage']) + 5, 1))),
#         ) +     
#         alt.Chart(df).mark_area(opacity=0.25, color=col).encode(
#             alt.X("Total_Gb"),
#             alt.Y("Lower"),
#             alt.Y2("Upper"),
#         ) 
#     )
    
#     return(chart)

# plt = ( complexity_curves(preseq, "cornflowerblue") + #Sample Complexity 
#         complexity_curves(benchmark, "grey") + #Benchmark Complexity
#         alt.Chart(breadthdf).mark_line(opacity=0.5, strokeDash=[5,5], color="black").encode( #Observed Breadth Line
#             alt.X("Total_Gb"),
#             alt.Y("Breadth")
#             ) +
#         alt.Chart(breadthdf.query("Total_Gb == 1.2")).mark_text( #Observed Breadth Text
#             text="Observed Breadth of Coverage: " + str(round(genfrac * 100, 1)) + "%", 
#             dy=15, size=18, align="left").encode(
#                 alt.X("Total_Gb"),
#                 alt.Y("Breadth"),
#                 ) + 
#         alt.Chart().mark_text(text=["Complexity curve and 95% CI (blue).", #Caption
#                                     "This library requires " + str(breadth5per) + "Gb sequencing for 5% breadth of coverage.",
#                                     "The grey represent the mean and standard deviation of the HGSVC benchmark dataset",
#                                     "(200 arbitrarily chosen libraries from Chaisson et al. 2019)"], 
#                                 x=0, y=dimension, dy=75, align='left', fontSize=18)
#         ).properties(
#             title=alt.TitleParams(libId, fontSize=20), height=dimension, width=dimension
#             ).configure_axis(
#                 labelFontSize=15, titleFontSize=20
#                 )