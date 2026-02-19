#!/usr/bin/env python

## PURPOSE: Export SCE positions in excel matrix.
## USAGE:   python 06c-bpr-sce.py <breakpointSummary>
## OUTPUT:  sce_breakpoints.xlsx

import pandas as pd

bpr_table = pd.read_table(sys.argv[1], sep=" ").query('genoT != "ww-cc" & genoT != "cc-ww"')
bpr_table['sample'] = bpr_table['filenames'].str.split('.processed').str[0]

samples_list = sorted(list(set(bpr_table['sample'])))
chr_list = [f"chr{i}" for i in range(1, 23)]

bpr_matrix = pd.DataFrame(index=samples_list, columns=chr_list)

for s in samples_list:
    for c in chr_list:
        q = bpr_table.query('seqnames == @c & sample == @s').reset_index(drop=True)
        
        if q.empty:
            continue
        elif len(q) == 1:
            row = q.iloc[0]
            bpr_matrix.loc[s,c] = f"{row['start']}-{row['end']}:{row['genoT']}"
        elif len(q) % 2 == 0: #even
            # print(s + ' ' + c + " even:" + str(len(q)))
            sce = []
            for i in range(int(len(q)/2)):
                row1 = q.iloc[i*2]
                row2 = q.iloc[i*2+1]
                size = row2['end'] - row1['start']
                if size >= 10000000: #10 million bases
                    sce.append(f"{row1['start']}-{row2['start']}:{row1['genoT']}")
                    
            if len(sce) != 0:
                bpr_matrix.loc[s,c] = ', '.join(sce)
        else: #odd
            diff_df = pd.DataFrame(index=range(int(len(q)/2)), columns=["size1", "sce1", "size2", "sce2"])
            for i in range(int(len(q)/2)):
                row0 = q.iloc[i*2]
                row1 = q.iloc[i*2+1]
                row2 = q.iloc[i*2+2]
                
                size1 = row1['end'] - row0['start']
                diff_df.loc[i, "size1"] = size1
                if size1 >= 10000000:
                    diff_df.loc[i, "sce1"] = f"{row0['start']}-{row1['end']}:{row0['genoT']}"
                
                size2 = row2['end'] - row1['start']
                diff_df.loc[i, "size2"] = size2
                if size2 >= 10000000:
                    diff_df.loc[i, "sce2"] =  f"{row1['start']}-{row2['end']}:{row1['genoT']}"

            if diff_df[['size1', 'size2']].sum(axis=0).idxmin() == "size1":
                row = q.iloc[q.index.max()]
                sce = list(diff_df["sce1"].dropna()) + [f"{row['start']}-{row['end']}:{row['genoT']}"]
                bpr_matrix.loc[s,c] = ', '.join(sce)
            else:
                row = q.iloc[q.index.min()]
                sce = list(diff_df["sce2"].dropna()) + [f"{row['start']}-{row['end']}:{row['genoT']}"]
                bpr_matrix.loc[s,c] = ', '.join(sce)
                
bpr_matrix.reset_index().to_excel("sce_breakpoints.xlsx", index=False)