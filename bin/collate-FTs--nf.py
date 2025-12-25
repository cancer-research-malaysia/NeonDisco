#!/usr/bin/env python3

import os
import re
import sys
import argparse
import polars as pl
from pathlib import Path
from typing import List, Tuple, Dict, Optional

def create_empty_output(sample_id: str, output_filename: str) -> None:
    """
    Create empty output files with the expected schema structure.
    
    Args:
        sample_id: The sample identifier
        output_filename: Output file prefix (no file extension)
    """
    # Define the expected schema
    empty_data = {
        "fusionTranscriptID": [],
        "fusionGenePair": [],
        "breakpointID": [],
        "5pStrand": [],
        "3pStrand": [],
        "originalTool": [],
        "sampleID": [],
        "sampleNum": [],
        "sampleNum_Padded": [],
        # Arriba columns
        "predictedEffect_ARR": [],
        "mutationType_ARR": [],
        "confidenceLabel_ARR": [],
        "readingFrame_ARR": [],
        "splitReadsTotal_ARR": [],
        "discordantReadPairs_ARR": [],
        "filteredReads_ARR": [],
        "peptideSequence_ARR": [],
        # FusionCatcher columns
        "predictedEffect_FC": [],
        "fusionPairAnnotation_FC": [],
        "splitReadsTotal_FC": [],
        "discordantReadPairs_FC": [],
        "longestAnchor_FC": [],
        # STARFusion columns
        "breakpointSpliceType_SF": [],
        "fusionPairAnnotation_SF": [],
        "splitReadsTotal_SF": [],
        "discordantReadPairs_SF": [],
        "largeAnchorSupport_SF": [],
        "FFPM_SF": []
    }
    
    # Create empty DataFrame without type casting to avoid schema conflicts
    empty_df = pl.DataFrame(empty_data)
    
    try:
        # Save as parquet and tsv
        empty_df.write_parquet(f"{output_filename}.parquet")
        empty_df.write_csv(f"{output_filename}.tsv", separator="\t", null_value='NA')
        
        print("Empty output files created successfully:")
        print(f"  - {output_filename}.parquet")
        print(f"  - {output_filename}.tsv")
        print("Pipeline can continue with downstream processing.")
        
    except Exception as e:
        print(f"Warning: Could not create empty output files: {str(e)}")
        print("Pipeline may encounter issues in downstream steps.")

def extract_sample_num(filename: str, tool_suffix: str) -> Optional[str]:
    pattern = rf'^(\d+)[TN]_{re.escape(tool_suffix)}\.tsv$'
    match = re.search(pattern, os.path.basename(filename))
    return match.group(1) if match else None


def check_file_validity(file_path: str, tool_name: str) -> Tuple[bool, str]:
    """
    Check if a file is valid for processing.
    
    Returns:
        Tuple[bool, str]: (is_valid, reason_if_invalid)
    """
    if not Path(file_path).exists():
        return False, f"File does not exist: {file_path}"
    
    try:
        # Check if file is readable and has content
        file_size = Path(file_path).stat().st_size
        if file_size == 0:
            return False, f"File is completely empty: {file_path}"
        
        # Try to read just the header to check file format
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) == 0:
            return False, f"File has no lines: {file_path}"
        
        if len(lines) == 1:
            return False, f"File only contains header (no data rows): {file_path}"
        
        # Check if we have the minimum expected columns for each tool
        header = lines[0].strip().split('\t')
        
        required_columns = {
            'Arriba': ['#gene1', 'gene2', 'breakpoint1', 'breakpoint2'],
            'FusionCatcher': ['Gene_1_symbol(5end_fusion_partner)', 'Gene_2_symbol(3end_fusion_partner)'],
            'STARFusion': ['LeftGene', 'RightGene', 'LeftBreakpoint', 'RightBreakpoint']
        }
        
        if tool_name in required_columns:
            missing_cols = [col for col in required_columns[tool_name] if col not in header]
            if missing_cols:
                return False, f"Missing required columns for {tool_name}: {missing_cols}"
        
        return True, ""
        
    except Exception as e:
        return False, f"Error reading file {file_path}: {str(e)}"

def wrangle_df(file_path: str, sample_id: str, sample_num: str, tool_name: str) -> Optional[pl.LazyFrame]:
    """
    Process input file based on the tool name and return a standardized lazy DataFrame.
    Returns None if processing fails.
    """
    try:
        lazy_df = pl.scan_csv(file_path, separator="\t", null_values=["NA", ".", ""]).fill_null("NA")   
        # Check if the lazy dataframe would be empty after collection
        # We do a quick shape check by collecting just the first row
        try:
            sample_row = lazy_df.limit(1).collect()
            if sample_row.height == 0:
                print(f"Warning: {tool_name} file {file_path} has no data rows after parsing.")
                return None
        except Exception as e:
            print(f"Warning: Could not validate data in {tool_name} file {file_path}: {str(e)}")
            return None
        
        # Shared base columns (Common to all tools)
        base_cols = [
            pl.lit(tool_name).alias("originalTool"),
            pl.lit(sample_id).alias("sampleID"),
            pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
            pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
        ]

        match tool_name:
            case 'Arriba':
                tool_specific = [
                    (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
                    (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
                    (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
                    (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("5pStrand"),
                    (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("3pStrand"),
                    (pl.col('site1') + "__" + pl.col('site2')).alias("predictedEffect_ARR"),
                    pl.col('type').alias("mutationType_ARR"),
                    pl.col('confidence').alias("confidenceLabel_ARR"),
                    pl.col('reading_frame').alias("reading_frame_ARR").alias("readingFrame_ARR"),
                    (pl.col('split_reads1').cast(pl.Int64) + pl.col('split_reads2').cast(pl.Int64)).cast(pl.Utf8).alias("splitReadsTotal_ARR"),
                    pl.col('discordant_mates').cast(pl.Utf8).alias("discordantReadPairs_ARR"),
                    pl.col('filters').alias("filteredReads_ARR"),
                    pl.col('peptide_sequence').alias("peptideSequence_ARR"),
                ]
                return lazy_df.select(base_cols + tool_specific)

            case 'FusionCatcher':
                g1 = pl.when(pl.col('Gene_1_symbol(5end_fusion_partner)') == "NA").then(pl.col('Gene_1_id(5end_fusion_partner)')).otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
                g2 = pl.when(pl.col('Gene_2_symbol(3end_fusion_partner)') == "NA").then(pl.col('Gene_2_id(3end_fusion_partner)')).otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
                fc_pred_effect = (pl.col('Predicted_effect').alias('predictedEffect_FC') if 'Predicted_effect' in lazy_df.collect_schema().names() else pl.lit("NA").cast(pl.Utf8).alias('predictedEffect_FC'))

                tool_specific = [
                    (g1 + "::" + g2 + '__' + pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("fusionTranscriptID"),
                    (g1 + "::" + g2).alias("fusionGenePair"),
                    (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("breakpointID"),
                    pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2).alias("5pStrand"),
                    pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2).alias("3pStrand"),
                    fc_pred_effect,
                    pl.col('Fusion_description').alias("fusionPairAnnotation_FC"),
                    pl.col('Spanning_unique_reads').cast(pl.Utf8).alias("splitReadsTotal_FC"),
                    pl.col('Spanning_pairs').cast(pl.Utf8).alias("discordantReadPairs_FC"),
                    pl.col('Longest_anchor_found').cast(pl.Utf8).alias("longestAnchor_FC"),
                ]
                return lazy_df.select(base_cols + tool_specific)
            
            case 'STARFusion':
                g1 = pl.when(pl.col('LeftGene').str.split("^").list.get(0) != "").then(pl.col('LeftGene').str.split("^").list.get(0)).otherwise(pl.col('LeftGene').str.split("^").list.get(1).str.split(".").list.get(0))
                g2 = pl.when(pl.col('RightGene').str.split("^").list.get(0) != "").then(pl.col('RightGene').str.split("^").list.get(0)).otherwise(pl.col('RightGene').str.split("^").list.get(1).str.split(".").list.get(0))
                bp1 = pl.col('LeftBreakpoint').str.replace(r'^chr', '').str.split(':').list.slice(0, 2).list.join(':')
                bp2 = pl.col('RightBreakpoint').str.replace(r'^chr', '').str.split(':').list.slice(0, 2).list.join(':')
                
                tool_specific = [
                    (g1 + "::" + g2 + '__' + bp1 + "-" + bp2).alias("fusionTranscriptID"),
                    (g1 + "::" + g2).alias("fusionGenePair"),
                    (bp1 + "-" + bp2).alias("breakpointID"),
                    pl.col('LeftBreakpoint').str.split(':').list.get(2).alias("5pStrand"),
                    pl.col('RightBreakpoint').str.split(':').list.get(2).alias("3pStrand"),
                    pl.col('SpliceType').alias("breakpointSpliceType_SF"),
                    pl.col('annots').alias("fusionPairAnnotation_SF"),
                    pl.col('JunctionReadCount').cast(pl.Utf8).alias("splitReadsTotal_SF"),
                    pl.col('SpanningFragCount').cast(pl.Utf8).alias("discordantReadPairs_SF"),
                    pl.col('LargeAnchorSupport').alias("largeAnchorSupport_SF"),
                    pl.col('FFPM').cast(pl.Utf8).alias("FFPM_SF"),
                ]
                return lazy_df.select(base_cols + tool_specific)
            
            case _:
                raise ValueError(f"Unsupported tool name: {tool_name}")
            
        
        # match tool_name:
        #     case 'Arriba':
        #         tool_specific = [
        #             (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
        #             (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
        #             (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
        #             (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("5pStrand"),
        #             (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("3pStrand"),
        #             (pl.col('site1') + "__" + pl.col('site2')).alias("predictedEffect_ARR"),
        #             pl.col('type').alias("mutationType_ARR"),
        #             pl.col('confidence').alias("confidenceLabel_ARR"),
        #             pl.col('reading_frame').alias("readingFrame_ARR"),
        #             # ALL numeric counts cast to Utf8 to match placeholders
        #             (pl.col('split_reads1').cast(pl.Int64) + pl.col('split_reads2').cast(pl.Int64)).cast(pl.Utf8).alias("splitReadsTotal_ARR"),
        #             pl.col('discordant_mates').cast(pl.Utf8).alias("discordantReadPairs_ARR"),
        #             pl.col('filters').alias("filteredReads_ARR"),
        #             pl.col('peptide_sequence').alias("peptideSequence_ARR"),
        #             # Placeholders for other tools
        #             pl.lit("NA").cast(pl.Utf8).alias("predictedEffect_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("fusionPairAnnotation_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("longestAnchor_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("breakpointSpliceType_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("fusionPairAnnotation_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("largeAnchorSupport_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("FFPM_SF")
        #         ]
        #         return lazy_df.select(base_cols + tool_specific)

        #     case 'FusionCatcher':
        #         g1 = pl.when(pl.col('Gene_1_symbol(5end_fusion_partner)') == "NA").then(pl.col('Gene_1_id(5end_fusion_partner)')).otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
        #         g2 = pl.when(pl.col('Gene_2_symbol(3end_fusion_partner)') == "NA").then(pl.col('Gene_2_id(3end_fusion_partner)')).otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
                
        #         # Dynamic check for Predicted_effect column
        #         fc_pred_effect = (
        #             pl.col('Predicted_effect').alias('predictedEffect_FC') 
        #             if 'Predicted_effect' in lazy_df.collect_schema().names() 
        #             else pl.lit("NA").cast(pl.Utf8).alias('predictedEffect_FC')
        #         )

        #         tool_specific = [
        #             (g1 + "::" + g2 + '__' + pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("fusionTranscriptID"),
        #             (g1 + "::" + g2).alias("fusionGenePair"),
        #             (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("breakpointID"),
        #             pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2).alias("5pStrand"),
        #             pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2).alias("3pStrand"),
        #             fc_pred_effect,
        #             pl.col('Fusion_description').alias("fusionPairAnnotation_FC"),
        #             pl.col('Spanning_unique_reads').cast(pl.Utf8).alias("splitReadsTotal_FC"),
        #             pl.col('Spanning_pairs').cast(pl.Utf8).alias("discordantReadPairs_FC"),
        #             pl.col('Longest_anchor_found').cast(pl.Utf8).alias("longestAnchor_FC"),
        #             # Arriba/SF Placeholders
        #             pl.lit("NA").cast(pl.Utf8).alias("predictedEffect_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("mutationType_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("confidenceLabel_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("readingFrame_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("filteredReads_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("peptideSequence_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("breakpointSpliceType_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("fusionPairAnnotation_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("largeAnchorSupport_SF"),
        #             pl.lit("NA").cast(pl.Utf8).alias("FFPM_SF")
        #         ]
        #         return lazy_df.select(base_cols + tool_specific)
            
        #     case 'STARFusion':
        #         g1 = pl.when(pl.col('LeftGene').str.split("^").list.get(0) != "").then(pl.col('LeftGene').str.split("^").list.get(0)).otherwise(pl.col('LeftGene').str.split("^").list.get(1).str.split(".").list.get(0))
        #         g2 = pl.when(pl.col('RightGene').str.split("^").list.get(0) != "").then(pl.col('RightGene').str.split("^").list.get(0)).otherwise(pl.col('RightGene').str.split("^").list.get(1).str.split(".").list.get(0))
        #         bp1 = pl.col('LeftBreakpoint').str.replace(r'^chr', '').str.split(':').list.slice(0, 2).list.join(':')
        #         bp2 = pl.col('RightBreakpoint').str.replace(r'^chr', '').str.split(':').list.slice(0, 2).list.join(':')
                
        #         tool_specific = [
        #             (g1 + "::" + g2 + '__' + bp1 + "-" + bp2).alias("fusionTranscriptID"),
        #             (g1 + "::" + g2).alias("fusionGenePair"),
        #             (bp1 + "-" + bp2).alias("breakpointID"),
        #             pl.col('LeftBreakpoint').str.split(':').list.get(2).alias("5pStrand"),
        #             pl.col('RightBreakpoint').str.split(':').list.get(2).alias("3pStrand"),
        #             pl.col('SpliceType').alias("breakpointSpliceType_SF"),
        #             pl.col('annots').alias("fusionPairAnnotation_SF"),
        #             pl.col('JunctionReadCount').cast(pl.Utf8).alias("splitReadsTotal_SF"),
        #             pl.col('SpanningFragCount').cast(pl.Utf8).alias("discordantReadPairs_SF"),
        #             pl.col('LargeAnchorSupport').alias("largeAnchorSupport_SF"),
        #             pl.col('FFPM').cast(pl.Utf8).alias("FFPM_SF"),
        #             # Arriba/FC Placeholders
        #             pl.lit("NA").cast(pl.Utf8).alias("predictedEffect_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("mutationType_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("confidenceLabel_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("readingFrame_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("filteredReads_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("peptideSequence_ARR"),
        #             pl.lit("NA").cast(pl.Utf8).alias("predictedEffect_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("fusionPairAnnotation_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("splitReadsTotal_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("discordantReadPairs_FC"),
        #             pl.lit("NA").cast(pl.Utf8).alias("longestAnchor_FC")
        #         ]
        #         return lazy_df.select(base_cols + tool_specific)
            
        #     case _:
        #         raise ValueError(f"Unsupported tool name: {tool_name}")

    except Exception as e:
        print(f"Error processing {tool_name}: {str(e)}")
        return None

def collate_fusion_data(sample_id: str, output_filename: str, input_files: List[Tuple[str, str]]) -> None:
    """
    Collate fusion transcript data from multiple fusion detection tools.
    
    Args:
        sample_id: The sample identifier
        output_filename: Output file prefix (no file extension)
        input_files: List of tuples containing (file_path, tool_suffix)
    """

    # Map tool suffixes to full tool names
    tool_name_map = {
        'arr': 'Arriba', 
        'fc': 'FusionCatcher',
        'sf': 'STARFusion'
    }

    print(f"Sample ID: {sample_id}")
    print(f"Input files: {input_files}")
    print(f"Output filename: {output_filename}")

    # Initialize an empty list to store the lazy DataFrames
    lazy_dfs = []
    valid_tools = []
    skipped_files = []

    for input_path, suffix in input_files:
        print(f'Checking file: {input_path}')
        
        tool_name = tool_name_map.get(suffix)
        if not tool_name:
            print(f"Warning: Unknown tool suffix '{suffix}'. Skipping file {input_path}")
            skipped_files.append((input_path, f"Unknown tool suffix: {suffix}"))
            continue

        # Validate file before processing
        is_valid, reason = check_file_validity(input_path, tool_name)
        if not is_valid:
            print(f"Skipping {tool_name} file: {reason}")
            skipped_files.append((input_path, reason))
            continue

        # Extract sample number
        sample_num = extract_sample_num(input_path, suffix)
        if not sample_num:
            reason = f"Could not extract sample number from filename"
            print(f"Warning: {reason}: {input_path}")
            skipped_files.append((input_path, reason))
            continue
            
        print(f'Processing {tool_name} TSV file... (sample ID: {sample_id})')

        # Create a lazy dataframe for each file
        lazy_df = wrangle_df(input_path, sample_id, sample_num, tool_name)
        if lazy_df is not None:
            lazy_dfs.append(lazy_df)
            valid_tools.append(tool_name)
            print(f"Successfully processed {tool_name} data")
        else:
            skipped_files.append((input_path, f"Failed to process {tool_name} data"))
    
    # Report on processing results
    print(f"\n=== Processing Summary ===")
    print(f"Successfully processed tools: {valid_tools}")
    if skipped_files:
        print(f"Skipped files:")
        for file_path, reason in skipped_files:
            print(f"  - {file_path}: {reason}")
    
    if not lazy_dfs:
        print("\n=== No Fusion Events Detected ===")
        print("No valid fusion data found across all input files. This could be because:")
        print("- All fusion callers produced no results (header-only files)")
        print("- All input files have format issues")
        print("- No input files exist at the specified paths")
        print("- All input files are missing required columns")
        print("\nCreating empty output files with proper structure...")
        
        # Create an empty DataFrame with the expected schema
        create_empty_output(sample_id, output_filename)
        return

    print(f"\nConcatenating lazy DataFrames from {len(lazy_dfs)} valid tools...")

    # try:
    #     combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)
        
    #     # Explicitly define these variables so the linter sees them
    #     all_categoricals = [
    #         "fusionTranscriptID", "fusionGenePair", "breakpointID", "5pStrand", "3pStrand",
    #         "originalTool", "sampleID", "predictedEffect_ARR", "mutationType_ARR",
    #         "confidenceLabel_ARR", "readingFrame_ARR", "predictedEffect_FC", 
    #         "fusionPairAnnotation_FC", "breakpointSpliceType_SF", "fusionPairAnnotation_SF"
    #     ]
        
    #     ints = ["sampleNum"]

    #     # Final collection with explicit variables
    #     results = combined_lazy_df.with_columns(
    #         [pl.col(col).cast(pl.Categorical) for col in all_categoricals if col in combined_lazy_df.collect_schema().names()]
    #     ).with_columns(
    #         [pl.col(col).cast(pl.Int64) for col in ints if col in combined_lazy_df.collect_schema().names()]
    #     ).collect()
            
    #     results.write_parquet(f"{output_filename}.parquet")
    #     results.write_csv(f"{output_filename}.tsv", separator="\t")
        
    # except Exception as e:
    #     print(f"Error during final processing: {str(e)}")
    #     create_empty_output(sample_id, output_filename)

    try:
        # FIX: Diagonal concatenation aligns columns by name regardless of order
        combined_lazy_df = pl.concat(lazy_dfs, how="diagonal", rechunk=True)
        
        all_categoricals = [
            "fusionTranscriptID", "fusionGenePair", "breakpointID", "5pStrand", "3pStrand",
            "originalTool", "sampleID", "predictedEffect_ARR", "mutationType_ARR",
            "confidenceLabel_ARR", "readingFrame_ARR", "predictedEffect_FC", 
            "fusionPairAnnotation_FC", "breakpointSpliceType_SF", "fusionPairAnnotation_SF"
        ]
        
        # Enforce consistent final column schema
        expected_cols = [
            "fusionTranscriptID", "fusionGenePair", "breakpointID", "5pStrand", "3pStrand", "originalTool", 
            "sampleID", "sampleNum", "sampleNum_Padded", "predictedEffect_ARR", "mutationType_ARR", 
            "confidenceLabel_ARR", "readingFrame_ARR", "splitReadsTotal_ARR", "discordantReadPairs_ARR", 
            "filteredReads_ARR", "peptideSequence_ARR", "predictedEffect_FC", "fusionPairAnnotation_FC", 
            "splitReadsTotal_FC", "discordantReadPairs_FC", "longestAnchor_FC", "breakpointSpliceType_SF", 
            "fusionPairAnnotation_SF", "splitReadsTotal_SF", "discordantReadPairs_SF", "largeAnchorSupport_SF", "FFPM_SF"
        ]

        schema_names = combined_lazy_df.collect_schema().names()
        missing_selects = [pl.lit("NA").alias(col) for col in expected_cols if col not in schema_names]

        results = (
            combined_lazy_df
            .with_columns(missing_selects)
            .select(expected_cols)
            .with_columns([pl.col(col).cast(pl.Categorical) for col in all_categoricals if col in expected_cols])
            .collect()
        )
            
        results.write_parquet(f"{output_filename}.parquet")
        results.write_csv(f"{output_filename}.tsv", separator="\t", null_value='NA')
        
    except Exception as e:
        print(f"Error during final processing: {str(e)}")
        create_empty_output(sample_id, output_filename)

def parse_arguments():
    """Parse command line arguments using argparse."""
    parser = argparse.ArgumentParser(
        description="Collate fusion transcript data from multiple fusion detection tools."
    )
    
    parser.add_argument("--sample_id", "-s", required=True, help="Sample identifier")
    parser.add_argument("--output", "-o", required=True, help="Output filename path (without extension)")
    
    # Add input file arguments
    input_group = parser.add_argument_group('input files')
    input_group.add_argument("--arriba", help="Path to Arriba output file")
    input_group.add_argument("--fusioncatcher", help="Path to FusionCatcher output file")
    input_group.add_argument("--starfusion", help="Path to STARFusion output file")
    
    # Add a way to provide additional inputs in the original positional format
    parser.add_argument("--inputs", nargs="+", help="Additional input files in format: file1 suffix1 file2 suffix2...")
    return parser.parse_args()

def main():
    """Main function to parse arguments and call the collate_fusion_data function."""
    # Support both new argparse style and legacy positional arguments
    if len(sys.argv) > 1 and sys.argv[1].startswith('-'):
        # Use argparse for modern argument parsing
        args = parse_arguments()
        
        # Prepare input files list
        input_files = []
        
        # Add specified tool files
        if args.arriba:
            print(f"Adding Arriba file: {args.arriba}")
            input_files.append((os.path.abspath(args.arriba), 'arr'))
        if args.fusioncatcher:
            print(f"Adding FusionCatcher file: {args.fusioncatcher}")
            input_files.append((os.path.abspath(args.fusioncatcher), 'fc'))
        if args.starfusion:
            print(f"Adding STARFusion file: {args.starfusion}")
            input_files.append((os.path.abspath(args.starfusion), 'sf'))
            
        # Add additional inputs if provided
        if args.inputs:
            if len(args.inputs) % 2 != 0:
                print("Error: Additional inputs must be provided as pairs of file path and tool suffix.")
                return  # Return instead of sys.exit(1)
                
            for i in range(0, len(args.inputs), 2):
                input_files.append((os.path.abspath(args.inputs[i]), args.inputs[i+1]))
        
        print(f"Input files: {input_files}")
        # Call the collate function
        try:
            collate_fusion_data(
                sample_id=args.sample_id,
                output_filename=args.output,
                input_files=input_files
            )
        except Exception as e:
            print(f"Unexpected error: {str(e)}")
            print("Attempting to create empty output files...")
            try:
                create_empty_output(args.sample_id, args.output)
            except Exception as e2:
                print(f"Failed to create empty output files: {str(e2)}")
                print("Pipeline may encounter issues in downstream steps.")

if __name__ == "__main__":
    main()
	