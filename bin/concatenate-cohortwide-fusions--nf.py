#!/usr/bin/env python3
"""
Script to concatenate multiple fusion TSV files into a cohort-wide file
"""
import sys
import argparse
from pathlib import Path
import polars as pl

def concatenate_cohortwide_fusions(input_files, output_file):
    """
    Concatenate multiple TSV files into a cohort-wide fusion file
    """
    try:
        # List to store all dataframes
        dfs = []

        # Process each input file
        for file_path in input_files:
            file_path = Path(file_path)
            if not file_path.exists():
                print(f"Warning: File {file_path} does not exist, skipping...")
                continue

            print(f"Processing: {file_path}")
            
            # Read TSV file
            df = pl.read_csv(file_path, separator='\t', schema_overrides={'sampleNum_Padded': pl.Utf8} )
            
            if len(df) == 0:
                print(f"Warning: File {file_path} is empty, skipping...")
                continue
            
            dfs.append(df)
            print(f"  - Added {len(df)} rows from {file_path.name}")
        
        if not dfs:
            print("Error: No valid TSV files found!")
            sys.exit(1)
        
        # Concatenate all dataframes
        print("Concatenating all dataframes...")
        cohort_df = pl.concat(dfs, how='vertical')
        
        # Sort by sampleNum_Padded if available, otherwise by sampleID
        if 'sampleNum_Padded' in cohort_df.columns:
            cohort_df = cohort_df.sort(['sampleNum_Padded', 'fusionGenePair'])
        elif 'sampleID' in cohort_df.columns:
            cohort_df = cohort_df.sort(['sampleID', 'fusionGenePair'])
        
        # Write the concatenated file
        cohort_df.write_csv(output_file, separator='\t')
        print(f"Cohort-wide TSV written to: {output_file}")
        print(f"Total rows in cohort file: {len(cohort_df)}")
        
        # Print summary statistics
        if 'sampleID' in cohort_df.columns:
            total_samples = cohort_df['sampleID'].n_unique()
            total_fusions = len(cohort_df)
            print(f"Summary: {total_samples} samples, {total_fusions} total fusions")
        
        return True
        
    except Exception as e:
        print(f"Error processing files: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Concatenate fusion TSV files')
    parser.add_argument('--input_files', nargs='+', required=True,
                       help='Input TSV files to concatenate')
    parser.add_argument('--output', required=True,
                       help='Output cohort TSV file')
    
    args = parser.parse_args()
    
    # Concatenate the files
    success = concatenate_cohortwide_fusions(args.input_files, args.output)
    
    if success:
        print("\nConcatenation completed successfully!")
        print(f"Processed {len(args.input_files)} input files")

if __name__ == "__main__":
    main()

