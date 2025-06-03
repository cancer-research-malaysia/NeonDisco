import os
import sys
import polars as pl

def main():
    """
    Process and filter fusion transcript data based on:
    1. Removing duplicate entries
    2. Adding tool overlap information
    3. Checking presence in CCLE and internal cell lines
    4. Filtering out fusions found in panel of normals
    5. Formatting for FusionInspector compatibility
    The script takes command-line arguments for input files and outputs the results in TSV format and a txt format for Fusion Inspector.
    """
    # Parse command-line arguments
    sample_name = sys.argv[1]
    parquet_input_file = os.path.abspath(sys.argv[2])
    panel_of_normals_file = os.path.abspath(sys.argv[3])
    ccle_internal_cell_line_file = os.path.abspath(sys.argv[4])
    output_filename = os.path.abspath(sys.argv[5])
    
    # Print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Parquet file of collated FTs: {parquet_input_file}")
    print(f"Parquet file of panel of Normals (TCGA Normals) FTs: {panel_of_normals_file}")
    print(f"Parquet file of CCLE + internal cell line FTs: {ccle_internal_cell_line_file}")
    print(f"Output filename: {output_filename}")

    # Load the collated fusion transcript data
    print("Loading collated fusion transcript data...")
    collated_df = pl.scan_parquet(parquet_input_file).collect()
    
    # Step 1: Filter for unique rows based on fusionTranscriptID and originalTool
    print("Filtering for unique rows based on fusionTranscriptID and originalTool...")
    unique_collated_df = collated_df.unique(subset=["fusionTranscriptID", "originalTool"])
    
    # Step 2: Create a new dataframe with unique fusionTranscriptIDs and list of tools that detected them
    print("Creating tool overlap information...")
    tool_group_df = (
        unique_collated_df
        .group_by('fusionTranscriptID')
        .agg(
            pl.col('originalTool').unique().alias('detectedBy'),
            pl.col('originalTool').unique().count().alias('toolOverlapCount')
        )
    )
    
    # Step 3: Join the unique fusions with the tool overlap information
    print("Joining unique fusions with tool overlap information...")
    unique_fusions_df = (
        unique_collated_df
        .drop('originalTool')
        .unique('fusionTranscriptID')
        .join(tool_group_df, on='fusionTranscriptID')
    )
    
    # Step 4: Load CCLE & Internal Cell Line FT data
    print("Loading CCLE & Internal Cell Line FT data...")
    ccle_df = pl.scan_parquet(ccle_internal_cell_line_file).collect()
    
    # Step 5: Add 'foundInCCLE&InternalCLs' column to unique fusions
    print("Adding 'foundInCCLE&InternalCLs' column to unique fusions...")
    ccle_set = set(ccle_df['breakpointID'].to_list())
    ccle_added_df = unique_fusions_df.with_columns(
        pl.when(pl.col('breakpointID').is_in(ccle_set)).then(True).otherwise(False).alias('foundInCCLE&InternalCLs')
    )
    
    # Step 6: Load Panel of Normals (TCGA Normals) data
    print("Loading Panel of Normals data...")
    pon_df = pl.scan_parquet(panel_of_normals_file).collect()
    
    # Step 7: Filter out breakpoints that appear in the Panel of Normals
    print("Filtering out breakpoints that appear in the Panel of Normals...")
    pon_set = set(pon_df['breakpointID'].to_list())
    normfilt_df = ccle_added_df.filter(~pl.col('breakpointID').is_in(pon_set))

    # Step 8: First cast the categorical column 'fusionGenePair' to string, then do the replacement
    print("Creating FusionInspector format column...")
    final_result_df = normfilt_df.with_columns(
        pl.col('fusionGenePair').cast(pl.Utf8).str.replace('::', '--').alias('fusionGenePair_FusIns')
    )
    
    # Save the result
    print(f"Saving filtered results to {output_filename}.tsv...")

    # Format detectedBy column to use " | " as separator between tools because polars represents the list as nested data
    export_df = final_result_df.with_columns([
        pl.col('detectedBy').list.eval(pl.element().cast(pl.Utf8)).list.join(" | ").alias('detectedBy')
    ])

    # write to tsv using polars
    export_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Results saved to {output_filename}.tsv")
    # write to txt file just the FusIns column
    export_df.select('fusionGenePair_FusIns').unique().write_csv(f"{output_filename}-unique-genePairs-for-FusIns.txt", include_header=False)
    print(f"Unique fusion gene pairs for FusionInspector saved to {output_filename}-unique-genePairs-for-FusIns.txt")
    
    print("Processing complete!")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: wrangle-and-filter-FTs--nf.py <sample id> <parquet file of combined FTs> <panel of normals FTs parquet file> <ccle+internal cell line FTs parquet file> <output filename>")
        sys.exit(1)
    main()
