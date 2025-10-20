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
    Optimized approach using separate dataframes per tool, then joining.
    This avoids row iteration while keeping the logic clear and correct.
    """
    # Get tool summary
    tool_summary = df.group_by(groupby_cols).agg([
        pl.col("originalTool").unique().sort().str.join(" | ").alias("detectedBy"),
        pl.col("originalTool").n_unique().alias("toolOverlapCount")
    ])
    
    # Separate by tool
    arriba_df = df.filter(pl.col("originalTool") == "Arriba")
    fc_df = df.filter(pl.col("originalTool") == "FusionCatcher")
    sf_df = df.filter(pl.col("originalTool") == "STAR-Fusion")
    
    # Get column lists from original column order
    value_cols = [col for col in original_columns if col not in groupby_cols and col != "originalTool"]
    arriba_cols = [col for col in value_cols if col.endswith("_ARR")]
    fc_cols = [col for col in value_cols if col.endswith("_FC")]
    sf_cols = [col for col in value_cols if col.endswith("_SF")]
    generic_cols = [col for col in value_cols if col not in arriba_cols + fc_cols + sf_cols]
    
    # Helper function to create safe filter for each column
    def safe_agg_expr(col):
        """Create aggregation expression that handles both string and numeric types."""
        col_type = df.select(pl.col(col)).dtypes[0]
        
        # For numeric columns, just take the first non-null value
        if col_type in [pl.Int8, pl.Int16, pl.Int32, pl.Int64, pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64, pl.Float32, pl.Float64]:
            return pl.col(col).drop_nulls().first().alias(col)
        
        # For string columns, filter out "NA" and nulls
        return (
            pl.when(pl.col(col).is_not_null())
            .then(
                pl.when(pl.col(col) != "NA")
                .then(pl.col(col))
                .otherwise(None)
            )
            .otherwise(None)
            .drop_nulls()
            .first()
            .alias(col)
        )
    
    # Helper function to fill nulls appropriately based on column type
    def fill_nulls_by_type(result_df, columns):
        """Fill nulls with 'NA' for string columns, keep null for numeric columns."""
        for col in columns:
            if col in result_df.columns:
                col_type = result_df.select(pl.col(col)).dtypes[0]
                if col_type == pl.Utf8:
                    result_df = result_df.with_columns(pl.col(col).fill_null("NA"))
                # For numeric columns, leave as null (will be written as empty in TSV)
        return result_df
    
    # Start with empty result containing only groupby columns
    result = tool_summary
    
    # Add Arriba columns
    if len(arriba_df) > 0 and len(arriba_cols) > 0:
        arriba_agg = arriba_df.group_by(groupby_cols).agg([
            safe_agg_expr(col) for col in arriba_cols
        ])
        result = result.join(arriba_agg, on=groupby_cols, how="left")
        result = fill_nulls_by_type(result, arriba_cols)
    else:
        for col in arriba_cols:
            result = result.with_columns(pl.lit("NA").alias(col))
    
    # Add FusionCatcher columns
    if len(fc_df) > 0 and len(fc_cols) > 0:
        fc_agg = fc_df.group_by(groupby_cols).agg([
            safe_agg_expr(col) for col in fc_cols
        ])
        result = result.join(fc_agg, on=groupby_cols, how="left")
        result = fill_nulls_by_type(result, fc_cols)
    else:
        for col in fc_cols:
            result = result.with_columns(pl.lit("NA").alias(col))
    
    # Add STAR-Fusion columns
    if len(sf_df) > 0 and len(sf_cols) > 0:
        sf_agg = sf_df.group_by(groupby_cols).agg([
            safe_agg_expr(col) for col in sf_cols
        ])
        result = result.join(sf_agg, on=groupby_cols, how="left")
        result = fill_nulls_by_type(result, sf_cols)
    else:
        for col in sf_cols:
            # Need to check original df to determine type
            if col in df.columns:
                col_type = df.select(pl.col(col)).dtypes[0]
                if col_type == pl.Utf8:
                    result = result.with_columns(pl.lit("NA").alias(col))
                else:
                    # For numeric columns, add as null
                    result = result.with_columns(pl.lit(None).cast(col_type).alias(col))
            else:
                result = result.with_columns(pl.lit("NA").alias(col))
    
    # Add generic columns (from any tool)
    if len(generic_cols) > 0:
        generic_agg = df.group_by(groupby_cols).agg([
            safe_agg_expr(col) for col in generic_cols
        ])
        result = result.join(generic_agg, on=groupby_cols, how="left")
        result = fill_nulls_by_type(result, generic_cols)
    
    # Reorder columns to match original order
    # Order: groupby_cols, detectedBy, toolOverlapCount, then all value_cols in original order
    final_column_order = groupby_cols + ['detectedBy', 'toolOverlapCount'] + value_cols
    result = result.select([col for col in final_column_order if col in result.columns])
    
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
    # Force sampleNum_Padded to be read as string to preserve zero padding
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
    