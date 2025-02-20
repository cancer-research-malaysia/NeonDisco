import os
import sys
import pandas as pd
import numpy as np

# Check if directory path is provided as argument
if len(sys.argv) != 2:
    print("Usage: python prune-SNAF.py <input directory>")
    sys.exit(1)

# Get base directory from command line argument
base_dir = sys.argv[1]

# construct paths
counts_path = os.path.join(base_dir, 'altanalyze_output/ExpressionInput/counts.original.txt')
output_full_path = os.path.join(base_dir, 'altanalyze_output/ExpressionInput/counts.original.full.txt')
event_annotation_path = os.path.join(base_dir, 'altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')
output_pruned_path = os.path.join(base_dir, 'altanalyze_output/ExpressionInput/counts.original.pruned.txt')

# preprocess the dataframe
df = pd.read_csv(counts_path, sep='\t', index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
df.to_csv(output_full_path, sep='\t')

# filter to EventAnnotation file
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv(event_annotation_path, sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]
df.to_csv(output_pruned_path, sep='\t')
