#!/usr/bin/env python3

import os
import re
import sys
import polars as pl
from pathlib import Path

def extract_sample_id(filename, suffix):
    # use the lines below if you want to keep the T in the sample ID
    # pattern = rf'^(\d+)(T)_{re.escape(suffix)}\.tsv$'
    # match = re.search(pattern, os.path.basename(filename))
    # return f"{match.group(1)}{match.group(2)}" if match else None

    pattern = rf'^(\d+)T_{re.escape(suffix)}\.tsv$'
    match = re.search(pattern, os.path.basename(filename))
    return match.group(1) if match else None

def extract_fuscat_breakpoint(s):
    return s.str.split(':').list.slice(0, 2).list.join(':')

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

def wrangle_df(file_path, sample_id, tool_name):
    lazy_df = pl.scan_csv(file_path, separator="\t")
    match tool_name:
        case 'Arriba':
            return lazy_df.select([
            (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
            (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
            (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
            (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("strand1"),
            (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("strand2"),
            'site1', 
            'site2', 
            'type', 
            'confidence',
            pl.lit(sample_id).alias("sampleID"),
            pl.lit(sample_id).cast(pl.Utf8).str.zfill(3).alias("sampleID_padded"),
            pl.lit(tool_name).alias("toolID")
            ])
        case 'FusionCatcher':
            # Handle NaN values in gene symbol columns by replacing with gene IDs
            gene1_expr = (
                pl.when(pl.col('Gene_1_symbol(5end_fusion_partner)').is_null() | (pl.col('Gene_1_symbol(5end_fusion_partner)') == ""))
                .then(pl.col('Gene_1_id(5end_fusion_partner)'))
                .otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
            )
            
            gene2_expr = (
                pl.when(pl.col('Gene_2_symbol(3end_fusion_partner)').is_null() | (pl.col('Gene_2_symbol(3end_fusion_partner)') == ""))
                .then(pl.col('Gene_2_id(3end_fusion_partner)'))
                .otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
            )
            
            base_columns = [
            (gene1_expr + "::" + gene2_expr + '__' + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("fusionTranscriptID"),
            (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
            (extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("breakpointID"),
            (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2)).alias("strand1"),
            (pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2)).alias("strand2")
            ]
            # Check if 'Predicted_effect' column exists
            if 'Predicted_effect' in lazy_df.collect_schema().names():
                # print(f'Sample ID: {sample_id} has Predicted_effect column')
                predicted_effect_columns = [
                    pl.col('Predicted_effect').str.extract(r'^([^/]+)(?:/|$)').alias('site1'),
                    pl.when(pl.col('Predicted_effect').str.contains('/'))
                        .then(pl.col('Predicted_effect').str.extract(r'/(.+)$'))
                        .otherwise(pl.lit('.'))
                        .alias('site2')
                    ]
            else:
                # print(f'Sample ID: {sample_id} does not have Predicted_effect column')
                predicted_effect_columns = [
                    pl.lit('.').alias("site1"),
                    pl.lit('.').alias("site2")
                    ]
            
            return lazy_df.select(base_columns + predicted_effect_columns + [pl.lit('.').alias("type"),
            pl.lit('.').alias("confidence"),
            pl.lit(sample_id).alias("sampleID"),
            pl.lit(sample_id).cast(pl.Utf8).str.zfill(3).alias("sampleID_padded"),
            pl.lit(tool_name).alias("toolID")]) 
        case _:
            raise ValueError(f"Unsupported tool name: {tool_name}")

def main():
    # this script would take these parameters:
    # combine-FTs-nf.py <sample id> <FT_arr.tsv file> <tool_suffix matching the file name> <FT_fc.tsv> <tool_suffix matching the file name>
    sample_name = sys.argv[1]
    input_tuples = [(os.path.abspath(sys.argv[2]), sys.argv[3]), (os.path.abspath(sys.argv[4]), sys.argv[5])]
    
    # print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Input arguments: {input_tuples}")

    # initialize an empty dictionary to store the lazy DataFrames
    lazy_dfs = []

    for input_path, suffix in input_tuples:
        if not Path(input_path).exists():
            print(f"Error: File {input_path} does not exist. Aborting...")
            sys.exit(1)
        
        print('Setting tool name...')
        tool_name = {'arr': 'Arriba', 'fc': 'FusionCatcher'}.get(suffix)

        # extract sample id
        sample_id = extract_sample_id(input_path, suffix)
        print(f'Reading {tool_name} of {sample_name} TSV file...(sample ID: {sample_id})')

        # Create a lazy dataframe for each file
        lazy_df = wrangle_df(input_path, sample_id, tool_name)

        # Append the lazy DataFrame to the list
        lazy_dfs.append(lazy_df)
    
    # Concatenate all lazy DataFrames
    print("Concatenating lazy DataFrames from Arriba and FusionCatcher...")
    combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)
    
    print("Concatenation completed. Collecting...")

    # Sort by Sample ID (padded), drop that column, then collect
    results = combined_lazy_df.sort("sampleID_padded").drop("sampleID_padded").with_columns(
        [pl.col(col).cast(pl.Categorical) for col in ['fusionTranscriptID', 'fusionGenePair', 'breakpointID', 'strand1', 'strand2', 'site1', 'site2', 'type', 'confidence', 'toolID']]).with_columns(
            pl.col("sampleID").cast(pl.Int64)).collect()

    # print(results)
    # save as parquet and tsv
    print(f"Saving as parquet and tsv files...")
    results.write_parquet(f"{sample_name}-combined-tool-FT-UNFILTERED.parquet")
    results.write_csv(f"{sample_name}-combined-tool-FT-UNFILTERED.tsv", separator="\t")

    print("Done.")

if __name__ == "__main__":
    main()

#############################################################
