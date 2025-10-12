#!/usr/bin/env python3

import sys
import polars as pl
from pathlib import Path

def postprocess_neopeptides(input_files, output_file):
    """
    Merge multiple TSV files with the same schema, standardize column names,
    add DeepImmuno column if missing, reorder columns, and sort.
    
    Args:
        input_files: List of input TSV file paths
        output_file: Output TSV file path
    """
    
    # Priority columns that should come first (in order)
    priority_cols = [
        "Sample_Name",
        "Gene_Name", 
        "HLA_Allele",
        "Sub-peptide_Position",
        "Epitope_Seq",
        "Variant_Type",
        "Median_IC50_Score",
        "Best_IC50_Score",
        "Best_IC50_Score_Method",
        "Median_Percentile",
        "Best_Percentile",
        "Best_Percentile_Method"
    ]
    
    dfs = []
    
    for file_path in input_files:
        print(f"Reading {file_path}...")
        
        try:
            df = pl.read_csv(file_path, separator='\t')
    
            if df.height == 0:
                print(f"  Warning: {file_path} has no data rows, skipping...")
                continue
        
        except Exception as e:
            print(f"  Warning: Could not read {file_path}: {e}, skipping...")
            continue
    
        print(f"  Columns: {len(df.columns)}, Rows: {df.height}")
        
        # Replace whitespace in column names with underscores
        df = df.rename({col: col.replace(' ', '_') for col in df.columns})
        
        # Fix the vexing " / " formatting in Chromosome, Start, and Stop columns
        for col in ['Chromosome', 'Start', 'Stop']:
            if col in df.columns:
                df = df.with_columns(
                    pl.col(col).str.replace_all(' / ', '::').alias(col)
                )
        
        # Add DeepImmuno column if it doesn't exist
        if "DeepImmuno_Score" not in df.columns:
            df = df.with_columns(pl.lit("NA").alias("DeepImmuno_Score"))
        
        dfs.append(df)
    
    if not dfs:
        print("Error: No valid dataframes to merge!")
        sys.exit(1)
    
    # Debug: Print schema information
    print(f"\nMerging {len(dfs)} dataframes...")
    for i, df in enumerate(dfs):
        print(f"  DataFrame {i}: {len(df.columns)} columns")
    
    # Check if all dataframes have the same columns
    all_columns = [set(df.columns) for df in dfs]
    if len(set(frozenset(cols) for cols in all_columns)) > 1:
        print("\nWarning: Column mismatch detected!")
        for i, cols in enumerate(all_columns):
            print(f"  DataFrame {i} columns: {len(cols)}")
        
        # Find columns that differ
        all_cols_set = set().union(*all_columns)
        for col in sorted(all_cols_set):
            present_in = [i for i, cols in enumerate(all_columns) if col in cols]
            if len(present_in) != len(dfs):
                print(f"  Column '{col}' missing in dataframes: {[i for i in range(len(dfs)) if i not in present_in]}")
    
    # Merge all dataframes
    try:
        merged_df = pl.concat(dfs, how="vertical_relaxed")
    except Exception as e:
        print(f"\nError during merge: {e}")
        print("Attempting alternative merge strategy...")
        # Align all dataframes to have the same columns
        all_cols = sorted(set().union(*[set(df.columns) for df in dfs]))
        aligned_dfs = []
        for df in dfs:
            for col in all_cols:
                if col not in df.columns:
                    df = df.with_columns(pl.lit("NA").alias(col))
            aligned_dfs.append(df.select(all_cols))
        merged_df = pl.concat(aligned_dfs, how="vertical")
    
    # Reorder columns: priority columns first, then the rest
    all_cols = merged_df.columns
    remaining_cols = [col for col in all_cols if col not in priority_cols]
    new_column_order = priority_cols + remaining_cols
    
    # Filter to only include columns that actually exist
    final_column_order = [col for col in new_column_order if col in all_cols]
    
    merged_df = merged_df.select(final_column_order)
    
    # Sort by Sample_Name (natural), then Gene_Name, then HLA_Allele (ascending)
    print("Sorting...")
    # Create natural sort key: extract numeric part and convert to int for proper sorting
    # Assumes format like "1T", "10T", etc.
    merged_df = merged_df.with_columns([
        pl.col("Sample_Name").str.extract(r'(\d+)', 1)
            .cast(pl.Int32, strict=False)
            .alias("_sample_num"),
        pl.col("Sample_Name").str.extract(r'(\D+)', 1)
            .alias("_sample_suffix")
    ])
    merged_df = merged_df.sort(["_sample_num", "_sample_suffix", "Gene_Name", "HLA_Allele"])
    merged_df = merged_df.drop(["_sample_num", "_sample_suffix"])
    
    # Write to output file
    print(f"Writing to {output_file}...")
    merged_df.write_csv(output_file, separator='\t')
    
    print(f"Done! Merged {len(dfs)} files with {len(merged_df)} total rows.")
    print(f"Output written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <output_file> <input_file1> <input_file2> ...")
        print("Example: python script.py merged_output.tsv file1.tsv file2.tsv file3.tsv")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    # Validate input files exist
    for f in input_files:
        if not Path(f).exists():
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    postprocess_neopeptides(input_files, output_file)
    