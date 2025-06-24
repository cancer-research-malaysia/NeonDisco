#!/usr/bin/env python3

import polars as pl
import sys
import os

def concatenate_parquet_files(file_paths, output_file="custom-PanelofNormals.parquet"):
    """
    Concatenate multiple parquet files with standardized schema using lazy evaluation.
    
    Args:
        file_paths (list): List of parquet file paths to concatenate
        output_file (str): name for the output parquet file
    """
    
    if not file_paths:
        print("No parquet files provided")
        sys.exit(1)
    
    # Filter out any non-parquet files and verify files exist
    valid_files = []
    for file_path in file_paths:
        if os.path.exists(file_path) and file_path.endswith('.parquet'):
            valid_files.append(file_path)
        else:
            print(f"Warning: Skipping {file_path} - file not found or not a parquet file")
    
    if not valid_files:
        print("No valid parquet files found")
        sys.exit(1)
    
    print(f"Found {len(valid_files)} parquet files to concatenate:")
    for f in sorted(valid_files):
        print(f"  - {f}")
    
    try:
        print("Processing files with lazy evaluation...")
        
        # Use lazy frame for memory efficiency
        lazy_df = pl.scan_parquet(sorted(valid_files))
        
        # Write parquet file directly without collecting (streaming)
        lazy_df.sink_parquet(output_file)
        print(f"Successfully created {output_file}")
        
        # Create TSV output file name
        tsv_file = output_file.replace('.parquet', '.tsv')
        
        # Write TSV file
        lazy_df.sink_csv(tsv_file, separator='\t')
        print(f"Successfully created {tsv_file}")
        
        # Get final stats efficiently
        schema = lazy_df.collect_schema()
        row_count = lazy_df.select(pl.len()).collect().item()
        
        print(f"Final files contain {row_count} rows and {len(schema)} columns")
        
    except Exception as e:
        print(f"Error during concatenation: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    # Get file paths from command line arguments
    if len(sys.argv) < 3:
        print("Usage: python3 concatenate_parquet_files.py <output_file> <file1.parquet> [file2.parquet] ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    file_paths = sys.argv[2:]  # This is a list slice, not a string
    
    concatenate_parquet_files(file_paths, output_file)

    