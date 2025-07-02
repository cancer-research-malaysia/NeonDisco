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
        'largeAnchorSupport_SF': [],
        'foundInCCLE&InternalCLs': [],
        'GenePairforFusInspector': []
    })
    
    # Write empty TSV
    empty_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Empty results file saved to {output_filename}.tsv")
    
    # Write empty txt file for FusionInspector
    with open(f"{output_filename}-unique-genePairs-for-FusInspector.txt", 'w') as f:
        pass  # Create empty file
    print(f"Empty fusion gene pairs file saved to {output_filename}-unique-genePairs-for-FusInspector.txt")

def main():
    """
    Process and filter fusion transcript data based on:
    1. Consolidating duplicate entries by:-
        - filling NA values across rows of tool-specific columns having similar fusion IDs)
        - removing duplicate entries
        - adding tool overlap information
    2. Filtering out fusions found in Panel of Normals (TCGA Normals or customized PoNs)
    3. Filtering out fusions found in Babiceanu et al. fusion file
    4. Adding a column indicating presence in CCLE and internal cell lines
    5. Adding a column indicating presence in Gao et al. TCGA recurrent fusion file
    6. Adding a column indicating presence in Mitelman et al. cancer gene fusion file
    7. Adding a column indicating presence in Klijn et al. cancer cell line fusion file
    8. Formatting for FusionInspector compatibility
    The script takes command-line arguments for input files and outputs the results in TSV format and a txt format for Fusion Inspector.
    """
    # Parse command-line arguments
    sample_name = sys.argv[1]
    parquet_input_file = os.path.abspath(sys.argv[2])
    panel_of_normals_file = os.path.abspath(sys.argv[3])
    babiceanu_et_al_normal_fusion_file = os.path.abspath(sys.argv[4])
    ccle_internal_cell_line_file = os.path.abspath(sys.argv[5])
    gao_et_al_fusion_file = os.path.abspath(sys.argv[6])
    mitelman_et_al_fusion_file = os.path.abspath(sys.argv[7])
    klijn_et_al_fusion_file = os.path.abspath(sys.argv[8])
    output_filename = os.path.abspath(sys.argv[9])
    
    # Print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Parquet file of collated FTs: {parquet_input_file}")
    print(f"Parquet file of panel of Normals (TCGA Normals) FTs: {panel_of_normals_file}")
    print(f"Parquet file of Babiceanu et al. normal tissue fusions: {babiceanu_et_al_normal_fusion_file}")
    print(f"Parquet file of CCLE + internal cell line FTs: {ccle_internal_cell_line_file}")
    print(f"Parquet file of Gao et al. fusions: {gao_et_al_fusion_file}")
    print(f"Parquet file of Mitelman et al. fusions: {mitelman_et_al_fusion_file}")
    print(f"Parquet file of Klijn et al. fusions: {klijn_et_al_fusion_file}")
    print(f"Output filename: {output_filename}")

    # Check if main input file is empty - if so, create empty outputs and exit
    if check_file_empty(parquet_input_file, "Collated fusion transcript data"):
        print(f"No fusion transcripts detected for sample {sample_name}. Creating empty output files...")
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

    # now sort the consolidated DataFrame by 'fusionTranscriptID'
    consolidated_df = consolidated_df.sort('fusionTranscriptID')

####################

    # Step 2: Load Panel of Normals (TCGA Normals) data
    print("Loading Panel of Normals data...")

    pon_df = pl.scan_parquet(panel_of_normals_file).collect()
    # pon_df.write_csv(f"{output_filename}-panel-of-normals.tsv", separator='\t', include_header=True)
    pon_set = set(pon_df['breakpointID'].to_list())

    # Filter out breakpoints that appear in the Panel of Normals
    print("Filtering out breakpoints that appear in the Panel of Normals...")
    normfilt_df = consolidated_df.filter(~pl.col('breakpointID').is_in(pon_set))

    # Check if filtering removed all fusions
    if normfilt_df.height == 0:
        print("WARNING: All fusions were filtered out by Panel of Normals filtering.")
        print("Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - all fusions filtered out!")
        return

    print(f"Panel of Normals filtering complete: {len(consolidated_df)} -> {len(normfilt_df)} rows")

################

    # Step 3: Load Babiceanu et al. normal tissue fusions
    print("Loading Babiceanu et al. normal tissue fusions...")
    babiceanu_df = pl.scan_parquet(babiceanu_et_al_normal_fusion_file).collect()

    babi_set = set(babiceanu_df['breakpointID'].to_list())

    # Filter out breakpoints that appear in Babiceanu et al. normal tissue fusions
    print("Filtering out breakpoints that appear in Babiceanu et al. normal tissue fusions...")
    normfilt_babifilt_df = normfilt_df.filter(~pl.col('breakpointID').is_in(babi_set))

    # Check if filtering removed all fusions
    if normfilt_babifilt_df.height == 0:
        print("WARNING: All fusions were filtered out by Babiceanu et al. normal tissue fusions filtering.")
        print("Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - all fusions filtered out!")
        return

    print(f"Babiceanu et al. filtering complete: {len(normfilt_df)} -> {len(normfilt_babifilt_df)} rows")

###################

    # Step 4: Load CCLE & Internal Cell Line FT data
    print("Loading CCLE & Internal Cell Line FT data...")

    ccle_df = pl.scan_parquet(ccle_internal_cell_line_file).collect()

    ccle_set = set(ccle_df['breakpointID'].to_list())

    # Add 'foundInCCLE&InternalCLs' column to unique fusions
    print("Adding 'foundInCCLE&InternalCLs' column to unique fusions...")
    ccle_added_df = normfilt_babifilt_df.with_columns(
        pl.when(pl.col('breakpointID').is_in(ccle_set)).then(True).otherwise(False).alias('foundInCCLE&InternalCLs')
    )

################
    # Step 5: Load Gao et al. TCGA recurrent fusion data
    print("Loading Gao et al. TCGA recurrent fusion data...")
    gao_df = pl.scan_parquet(gao_et_al_fusion_file).collect()

    gao_set = set(gao_df['breakpointID'].to_list())

    # Add 'foundInGaoetalTCGARecurrent' column to unique fusions
    print("Adding 'foundInGaoTCGARecurrent' column to unique fusions...")
    gao_added_df = ccle_added_df.with_columns(
        pl.when(pl.col('breakpointID').is_in(gao_set)).then(True).otherwise(False).alias('foundInGaoTCGARecurrent')
    )

################
    # Step 6: Load Mitelman et al. cancer gene fusion data
    print("Loading Mitelman et al. cancer gene fusion data...")
    mitelman_df = pl.scan_parquet(mitelman_et_al_fusion_file).collect()

    mitelman_set = set(mitelman_df['fusionGenePair'].to_list())

    # Add 'foundInMitelmanCancerFusions' column to unique fusions
    print("Adding 'foundInMitelmanCancerFusions' column to unique fusions...")
    mitelman_added_df = gao_added_df.with_columns(
        pl.when(pl.col('fusionGenePair').is_in(mitelman_set)).then(True).otherwise(False).alias('foundInMitelmanCancerFusions')
    )

################
    # Step 7: Load Klijn et al. cancer cell line fusion data
    print("Loading Klijn et al. cancer cell line fusion data...")
    klijn_df = pl.scan_parquet(klijn_et_al_fusion_file).collect()

    klijn_set = set(klijn_df['breakpointID'].to_list())

    # Add 'foundInKlijnCancerCellLines' column to unique fusions
    print("Adding 'foundInKlijnCancerCellLineFusions' column to unique fusions...")
    klijn_added_df = mitelman_added_df.with_columns(
        pl.when(pl.col('breakpointID').is_in(klijn_set)).then(True).otherwise(False).alias('foundInKlijnCancerCellLineFusions')
    )   

################

    # Step 8: Create FusionInspector format column
    print("Creating FusionInspector format column...")
    final_result_df = klijn_added_df.with_columns(
        pl.col('fusionGenePair').cast(pl.Utf8).str.replace('::', '--').alias('GenePairforFusInspector')
    )
    
##############

    # Step 9: Apply consensus filtering
    print("Applying consensus filtering...")

    # Filter for rows based on 'toolOverlapCount'
    # > 0 is default to keep the union of all fusions detected by at least one tool
    # change this value for different filtering
    export_consensus_df = final_result_df.filter(pl.col('toolOverlapCount') > 0)

    # Check if consensus filtering removed all fusions
    if export_consensus_df.height == 0:
        print("WARNING: No fusions passed the consensus filter.")
        print("Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - no consensus fusions found!")
        return
    else:
        # sort by toolOverlapCount and fusionGenePair
        export_consensus_df_sorted = export_consensus_df.sort(['toolOverlapCount', 'fusionGenePair'])

###########

    # Step 6: Save results
    print(f"Saving filtered results to {output_filename}.tsv...")
    export_consensus_df_sorted.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Results saved to {output_filename}.tsv")
    
    # Write unique fusion gene pairs for FusionInspector
    export_consensus_df_sorted.select('GenePairforFusInspector').unique().write_csv(f"{output_filename}-unique-genePairs-for-FusInspector.txt", include_header=False)
    print(f"Unique fusion gene pairs for FusionInspector saved to {output_filename}-unique-genePairs-for-FusInspector.txt")
    
    print("Processing complete!")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: wrangle-and-filter-FTs--nf.py <sample id> <parquet file of combined FTs> <panel of normals FTs parquet file> <Babiceanu et al normal fusions parquet file> <ccle+internal cell line FTs parquet file> <Gao et al fusion parquet file> <Mitelman et al fusion parquet file> <Klijn et al fusion parquet file> <output filename>")
        sys.exit(1)
    main()
