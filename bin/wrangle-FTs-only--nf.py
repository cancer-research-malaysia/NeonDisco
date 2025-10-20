#!/usr/bin/env python3

import os
import sys
import polars as pl

def check_file_empty(file_path, file_description):
    """
    Check if a parquet file is empty or contains no data rows.
    Returns True if empty, False if it contains data.
    """
    try:
        df = pl.scan_parquet(file_path).collect()
        if df.height == 0:
            print(f"WARNING: {file_description} is empty (no data rows)")
            return True
        return False
    except Exception as e:
        print(f"ERROR: Could not read {file_description}: {e}")
        return True

def merge_by_tool_suffixes(df, groupby_cols, original_columns):
    """
    Consolidate rows by fusion ID, combining tool-specific data.
    For each fusion, take the first non-null, non-"NA" value from any tool's row.
    """
    # Get tool summary
    tool_summary = df.group_by(groupby_cols).agg([
        pl.col("originalTool").unique().sort().str.join(" | ").alias("detectedBy"),
        pl.col("originalTool").n_unique().alias("toolOverlapCount")
    ])
    
    # For each column (except groupby and originalTool), aggregate by taking first non-null, non-"NA" value
    value_cols = [col for col in original_columns if col not in groupby_cols and col != "originalTool"]
    
    agg_exprs = []
    for col in value_cols:
        col_type = df.select(pl.col(col)).dtypes[0]
        
        if col_type in [pl.Int8, pl.Int16, pl.Int32, pl.Int64, 
                       pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64, 
                       pl.Float32, pl.Float64]:
            # For numeric columns: take first non-null
            agg_exprs.append(pl.col(col).drop_nulls().first().alias(col))
        else:
            # For string columns: take first that is not null and not "NA"
            agg_exprs.append(
                pl.col(col)
                .filter((pl.col(col).is_not_null()) & (pl.col(col) != "NA"))
                .first()
                .alias(col)
            )
    
    # Aggregate all value columns
    value_agg = df.group_by(groupby_cols).agg(agg_exprs)
    
    # Join with tool summary
    result = tool_summary.join(value_agg, on=groupby_cols, how="left")
    
    # Reorder columns: groupby cols, then detectedBy, toolOverlapCount, then value cols in original order
    sample_cols = ['sampleID', 'sampleNum', 'sampleNum_Padded']
    sample_cols = [col for col in sample_cols if col in value_cols]
    
    arriba_cols = [col for col in value_cols if col.endswith("_ARR")]
    fc_cols = [col for col in value_cols if col.endswith("_FC")]
    sf_cols = [col for col in value_cols if col.endswith("_SF")]
    other_cols = [col for col in value_cols 
                  if col not in sample_cols + arriba_cols + fc_cols + sf_cols]
    
    final_order = (groupby_cols + 
                   ['detectedBy', 'toolOverlapCount'] + 
                   sample_cols + 
                   arriba_cols + 
                   fc_cols + 
                   sf_cols + 
                   other_cols)
    
    result = result.select([col for col in final_order if col in result.columns])
    
    return result

def create_empty_output_files(output_filename):
    """
    Create empty output files with proper headers when no fusions are found.
    """
    empty_df = pl.DataFrame({
        'fusionTranscriptID': [],
        'fusionGenePair': [],
        'breakpointID': [],
        '5pStrand': [],
        '3pStrand': [],
        'detectedBy': [],
        'toolOverlapCount': [],
        'sampleID': [],
        'sampleNum': [],
        'sampleNum_Padded': [],
        '5pSite_ARR': [],
        '3pSite_ARR': [],
        'mutationType_ARR': [],
        'confidenceLabel_ARR': [],
        'predictedEffect_FC': [],
        'fusionPairAnnotation_FC': [],
        'fusionPairAnnotation_SF': [],
        'junctionReadCount_SF': [],
        'spanningFragCount_SF': [],
        'largeAnchorSupport_SF': []
    })
    
    # Write empty TSV
    empty_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Empty results file saved to {output_filename}.tsv")
    
def main():
    """
    Process fusion transcript data based on:
    1. Consolidating duplicate entries by:-
        - filling NA values across rows of tool-specific columns having similar fusion IDs)
        - removing duplicate entries
        - adding tool overlap information
    
    The script takes command-line arguments for input files and outputs the results in TSV format and a txt format for Fusion Inspector.
    """
    # Parse command-line arguments
    sample_name = sys.argv[1]
    parquet_input_file = os.path.abspath(sys.argv[2])
    output_filename = os.path.abspath(sys.argv[3])
    
    # Print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Parquet file of collated FTs: {parquet_input_file}")
    print(f"Output filename: {output_filename}")

    # Check if main input file is empty - if so, create empty outputs and exit
    if check_file_empty(parquet_input_file, "Collated fusion transcript data"):
        print(f"No fusion transcripts detected for sample {sample_name}. Creating empty output file...")
        create_empty_output_files(output_filename)
        print("Processing complete - no fusions detected!")
        return
    
    # Load the collated fusion transcript data
    print("Loading collated fusion transcript data...")
    collated_df = pl.read_parquet(parquet_input_file)
    
    # Ensure sampleNum_Padded is string type to preserve zero padding
    if 'sampleNum_Padded' in collated_df.columns:
        collated_df = collated_df.with_columns(
            pl.col('sampleNum_Padded').cast(pl.Utf8)
        )
    
    # Store original column order
    original_columns = collated_df.columns
    
    # Step 1: Consolidate duplicate rows and reshape the DataFrame
    print("Consolidating duplicate rows...")
    # Define the columns to group by
    groupby_cols = ['fusionTranscriptID', 'fusionGenePair', 'breakpointID', '5pStrand', '3pStrand']
    consolidated_df = merge_by_tool_suffixes(collated_df, groupby_cols, original_columns)
    
    if not consolidated_df.is_empty():
        print(f"Consolidation complete: {len(collated_df)} -> {len(consolidated_df)} rows")

    # sort by toolOverlapCount and fusionGenePair
    consolidated_df = consolidated_df.sort(['toolOverlapCount', 'fusionGenePair', 'sampleNum_Padded'], descending=[True, False, False])
    
    # Step 6: Save results with proper null handling
    print(f"Saving filtered results to {output_filename}.tsv...")
    # Use null_value parameter to write nulls as "NA" in the TSV for better compatibility
    consolidated_df.write_csv(f"{output_filename}.tsv", separator='\t', null_value='NA')
    print(f"Results saved to {output_filename}.tsv")
    
    print("Processing complete!")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: wrangle-FTs-only--nf.py <sample id> <parquet file of combined FTs> <output filename>")
        sys.exit(1)
    main()
    