#!/usr/bin/env python3

import os
import sys
import polars as pl
import pandas as pd

def main():
    # this script would take these parameters:
    # filter-ft-nf.py <sample id> <parquet file of combined FTs> <panel of normals tsv file>
    sample_name = sys.argv[1]
    parquet_input_file = os.path.abspath(sys.argv[2])
    panel_of_normals_file = os.path.abspath(sys.argv[3])
    
    # print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Parquet file of combined FTs: {parquet_input_file}")
    print(f"Location of panel-of-normal tsv file: {panel_of_normals_file}")

    #### 1. Get unique breakpointIDs from the combined ft dataframe based on toolID and add detection information
    # Load the combined parquet file
    combined_pq_df = pl.scan_parquet(parquet_input_file).collect().to_pandas(use_pyarrow_extension_array=True)
    
    # Split the dataframe by toolID
    arriba_df = combined_pq_df[combined_pq_df['toolID'] == 'Arriba'].copy()
    fusioncatcher_df = combined_pq_df[combined_pq_df['toolID'] == 'FusionCatcher'].copy()

    # Drop duplicates by breakpointID for each tool separately
    arriba_unique_df = arriba_df.drop_duplicates(subset=['breakpointID']).copy()
    fusioncatcher_unique_df = fusioncatcher_df.drop_duplicates(subset=['breakpointID']).copy()

    # Get sets of breakpointIDs from each tool
    arriba_breakpoints = set(arriba_unique_df['breakpointID'])
    fusioncatcher_breakpoints = set(fusioncatcher_unique_df['breakpointID'])

    # Find breakpoints detected by both tools
    common_breakpoints = arriba_breakpoints.intersection(fusioncatcher_breakpoints)

    # Add detection information using assign() with the new column name
    arriba_unique_df = arriba_unique_df.assign(detectedBy='Arriba')
    fusioncatcher_unique_df = fusioncatcher_unique_df.assign(detectedBy='FusionCatcher')

    # Mark breakpoints found by both tools
    mask = arriba_unique_df['breakpointID'].isin(common_breakpoints)
    arriba_unique_df.loc[mask, 'detectedBy'] = 'Both'

    # Only keep FusionCatcher rows that weren't detected by Arriba
    fusioncatcher_exclusive_df = fusioncatcher_unique_df[~fusioncatcher_unique_df['breakpointID'].isin(arriba_breakpoints)]

    # Combine Arriba rows with FusionCatcher-exclusive rows
    combined_unique_df = pd.concat([arriba_unique_df, fusioncatcher_exclusive_df])

    # 2. Filter out breakpoints seen in TCGA normals
    # load up the TCGA normals (all Arr and FC) dataframe
    panel_normals_uniq_df = pd.read_csv(panel_of_normals_file, sep='\t')

    # Get the unique breakpointIDs from the panel normals
    panel_normals_uniq_fts = set(panel_normals_uniq_df['breakpointID'])
    panel_normals_uniq_fts

    # Filter out breakpoints seen in TCGA normals
    final_unique_df = combined_unique_df[~combined_unique_df['breakpointID'].isin(panel_normals_uniq_fts)]

    # Save the result
    final_unique_df.to_csv(f"{sample_name}-combined-tool-FT-FILTERED.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()

#############################################################
