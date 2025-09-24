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

# Alternative approach to merge tool-specific information using explicit wrangling logic
def merge_by_tool_suffixes(df, groupby_cols):
    # Get unique fusions with their tools
    fusion_summary = df.group_by(groupby_cols).agg([
        pl.col("originalTool").unique().alias("detectedBy"),
        pl.col("originalTool").unique().count().alias("toolOverlapCount")
    ])

    result_rows = []

    for fusion_row in fusion_summary.iter_rows(named=True):
        # Start with fusion ID columns
        merged_row = {col: fusion_row[col] for col in groupby_cols}
        tools = fusion_row["detectedBy"]
        tool_count = fusion_row["toolOverlapCount"]
        merged_row["detectedBy"] = " | ".join(tools)
        merged_row["toolOverlapCount"] = tool_count

        # Get data for this fusion
        fusion_filter = pl.all_horizontal([
            pl.col(col) == fusion_row[col] for col in groupby_cols
        ])
        fusion_data = df.filter(fusion_filter)

        # For each column, get data from appropriate tool
        for col in df.columns:
            if col not in groupby_cols and col != "originalTool":
                # Check if this is a tool-specific column
                if col.endswith("_ARR") and "Arriba" in tools:
                    # Get Arriba data
                    arriba_data = fusion_data.filter(pl.col("originalTool") == "Arriba")
                    if len(arriba_data) > 0:
                        values = arriba_data.select(col).to_series().to_list()
                        non_null_values = [v for v in values if v is not None and v != "NA"]
                        merged_row[col] = non_null_values[0] if non_null_values else "NA"
                    else:
                        merged_row[col] = "NA"
                elif col.endswith("_FC") and "FusionCatcher" in tools:
                    # Get FusionCatcher data
                    fc_data = fusion_data.filter(pl.col("originalTool") == "FusionCatcher")
                    if len(fc_data) > 0:
                        values = fc_data.select(col).to_series().to_list()
                        non_null_values = [v for v in values if v is not None and v != "NA"]
                        merged_row[col] = non_null_values[0] if non_null_values else "NA"
                    else:
                        merged_row[col] = "NA"
                elif col.endswith("_SF") and "STAR-Fusion" in tools:
                    # Get STAR-Fusion data
                    sf_data = fusion_data.filter(pl.col("originalTool") == "STAR-Fusion")
                    if len(sf_data) > 0:
                        values = sf_data.select(col).to_series().to_list()
                        non_null_values = [v for v in values if v is not None and v != "NA"]
                        merged_row[col] = non_null_values[0] if non_null_values else "NA"
                    else:
                        merged_row[col] = "NA"
                else:
                    # For non-tool-specific columns, take first non-null
                    values = fusion_data.select(col).to_series().to_list()
                    non_null_values = [v for v in values if v is not None and v != "NA"]
                    merged_row[col] = non_null_values[0] if non_null_values else "NA"

        result_rows.append(merged_row)

    return pl.DataFrame(result_rows)

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
    collated_df = pl.scan_parquet(parquet_input_file).collect()
    
####################

    # Step 1: Consolidate duplicate rows and reshape the DataFrame
    print("Consolidating duplicate rows...")
    # Define the columns to group by
    groupby_cols = ['fusionTranscriptID', 'fusionGenePair', 'breakpointID', '5pStrand', '3pStrand']
    consolidated_df = merge_by_tool_suffixes(collated_df, groupby_cols)
    
    if not consolidated_df.is_empty():
        print(f"Consolidation complete: {len(collated_df)} -> {len(consolidated_df)} rows")

    # sort by toolOverlapCount and fusionGenePair
    consolidated_df = consolidated_df.sort(['toolOverlapCount', 'fusionGenePair', 'sampleNum_Padded'], descending=[True, False, False])
    
###########

    # Step 6: Save results
    print(f"Saving filtered results to {output_filename}.tsv...")
    consolidated_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Results saved to {output_filename}.tsv")
    
    print("Processing complete!")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: wrangle-FTs-only--nf.py <sample id> <parquet file of combined FTs> <output filename>")
        sys.exit(1)
    main()

