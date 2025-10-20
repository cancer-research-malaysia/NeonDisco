#!/usr/bin/env python3
"""
Script to concatenate multiple fusion TSV files into a cohort-wide file
Fixed version with better schema handling
"""
import sys
import argparse
from pathlib import Path
import polars as pl
from collections import defaultdict

def analyze_schemas(input_files):
    """
    Analyze schemas of all input files to determine consistent types
    """
    column_types = defaultdict(set)
    
    for file_path in input_files:
        file_path = Path(file_path)
        if not file_path.exists():
            continue
            
        try:
            # Read just the schema without loading all data
            # Use null_values to properly interpret "NA" as null
            df = pl.read_csv(
                file_path, 
                separator='\t', 
                n_rows=1,
                null_values=["NA", ""]
            )
            for col_name, dtype in df.schema.items():
                column_types[col_name].add(dtype)
        except Exception as e:
            print(f"Warning: Could not read schema from {file_path}: {e}")
            continue
    
    # Create schema overrides for problematic columns
    schema_overrides = {}
    
    # Force sampleNum_Padded to be string to preserve zero padding
    if 'sampleNum_Padded' in column_types:
        schema_overrides['sampleNum_Padded'] = pl.Utf8
    
    for col_name, types in column_types.items():
        if col_name in schema_overrides:
            continue  # Skip if already handled
            
        if len(types) > 1:
            # Handle type conflicts by choosing the most compatible type
            if pl.Utf8 in types:
                # If any file has string type, convert all to string
                schema_overrides[col_name] = pl.Utf8
            elif pl.Float64 in types and pl.Int64 in types:
                # If mixing int and float, use float
                schema_overrides[col_name] = pl.Float64
            else:
                # Default to string for safety
                schema_overrides[col_name] = pl.Utf8
    
    return schema_overrides

def concatenate_cohortwide_fusions(input_files, output_file):
    """
    Concatenate multiple TSV files into a cohort-wide fusion file
    """
    try:
        # First, analyze schemas to determine overrides
        print("Analyzing file schemas...")
        schema_overrides = analyze_schemas(input_files)
        
        if schema_overrides:
            print("Found schema conflicts or special columns, applying overrides:")
            for col, dtype in schema_overrides.items():
                print(f"  {col}: {dtype}")
            print()
        
        # List to store all dataframes
        dfs = []

        # Process each input file
        for file_path in input_files:
            file_path = Path(file_path)
            if not file_path.exists():
                print(f"Warning: File {file_path} does not exist, skipping...")
                continue

            print(f"Processing: {file_path}")
            
            try:
                # Read TSV file with schema overrides
                # Use null_values to interpret "NA" as null for proper type inference
                df = pl.read_csv(
                    file_path, 
                    separator='\t', 
                    schema_overrides=schema_overrides,
                    null_values=["NA", ""]
                )
                
                if len(df) == 0:
                    print(f"Warning: File {file_path} is empty, skipping...")
                    continue
                
                dfs.append(df)
                print(f"  - Added {len(df)} rows from {file_path.name}")
                
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                continue
        
        if not dfs:
            print("Error: No valid TSV files found!")
            sys.exit(1)
        
        # Concatenate all dataframes
        print("Concatenating all dataframes...")
        cohort_df = pl.concat(dfs, how='vertical_relaxed')

        # Sort by sampleNum_Padded if available, otherwise by sampleID
        sort_cols = []
        sort_desc = []
        
        if 'sampleNum_Padded' in cohort_df.columns:
            sort_cols.extend(['sampleNum_Padded', 'toolOverlapCount', 'fusionTranscriptID'])
            sort_desc.extend([False, True, False])
        elif 'sampleID' in cohort_df.columns:
            sort_cols.extend(['sampleID', 'toolOverlapCount', 'fusionTranscriptID'])
            sort_desc.extend([False, True, False])
        
        if sort_cols:
            cohort_df = cohort_df.sort(sort_cols, descending=sort_desc)
        
        if 'fusionTranscriptID' in cohort_df.columns:
            print("\nGenerating manifest...")
            
            # Extract unique IDs from the final dataframe
            fusion_ids = cohort_df['fusionTranscriptID'].unique().sort().to_list()
            
            # Determine the output filename for the manifest
            output_path = Path(output_file)
            base_name = output_path.stem 
            
            # Construct the new manifest filename: [basename]-UNIQUE.manifest.txt
            manifest_filename = f"{base_name}-UNIQUE.manifest.txt"
            manifest_file = output_path.with_name(manifest_filename)
            
            # Write the list of unique IDs to the file, one per row
            with open(manifest_file, 'w') as f:
                for id_val in fusion_ids:
                    f.write(f"{id_val}\n")
            
            print(f"Manifest written to: {manifest_file}")
            print(f"Total unique fusion IDs: {len(fusion_ids)}")
        else:
            print("Warning: 'fusionTranscriptID' column not found, skipping manifest generation.")

        print("\nWriting cohort-wide TSV...")
        # Write the concatenated file with NA for null values
        cohort_df.write_csv(output_file, separator='\t', null_value='NA')
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
        import traceback
        traceback.print_exc()
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
    