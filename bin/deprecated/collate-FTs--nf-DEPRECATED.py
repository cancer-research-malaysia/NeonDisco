#!/usr/bin/env python3

import os
import re
import sys
import argparse
import polars as pl
from pathlib import Path
from typing import List, Tuple, Dict, Optional

def extract_sample_num(filename: str, tool_suffix: str) -> Optional[str]:
    pattern = rf'^(\d+)[TN]_{re.escape(tool_suffix)}\.tsv$'
    match = re.search(pattern, os.path.basename(filename))
    return match.group(1) if match else None

def wrangle_df(file_path: str, sample_id: str, sample_num: str, tool_name: str) -> pl.LazyFrame:
    """Process input file based on the tool name and return a standardized lazy DataFrame."""
    # Read with null replacement
    lazy_df = pl.scan_csv(file_path, separator="\t").fill_null("NA")
    
    match tool_name:
        case 'Arriba':
            # Base columns
            base_columns = [
                (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
                (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
                (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
                (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("5pStrand"),
                (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("3pStrand"),
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ]
            
            # Arriba-specific columns
            arriba_columns = [
                (pl.col('site1')).alias("5pSite_ARR"),
                (pl.col('site2')).alias("3pSite_ARR"), 
                (pl.col('type')).alias("mutationType_ARR"),
                (pl.col('confidence')).alias("confidenceLabel_ARR")
            ]
            
            # Placeholder columns for other tools
            fc_columns = [
                pl.lit("NA").cast(pl.String).alias("predictedEffect_FC"),
                pl.lit("NA").cast(pl.String).alias("fusionPairAnnotation_FC")
            ]
            
            sf_columns = [
                pl.lit("NA").cast(pl.String).alias("fusionPairAnnotation_SF"),
                pl.lit("NA").cast(pl.String).alias("junctionReadCount_SF"),
                pl.lit("NA").cast(pl.String).alias("spanningFragCount_SF"),
                pl.lit("NA").cast(pl.String).alias("largeAnchorSupport_SF")
            ]
            
            return lazy_df.select(base_columns + arriba_columns + fc_columns + sf_columns)
            
        case 'FusionCatcher':
            # Handle NaN values in gene symbol columns by replacing with gene IDs
            gene1_expr = (
                pl.when((pl.col('Gene_1_symbol(5end_fusion_partner)') == "NA"))
                .then(pl.col('Gene_1_id(5end_fusion_partner)'))
                .otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
            )
            
            gene2_expr = (
                pl.when((pl.col('Gene_2_symbol(3end_fusion_partner)') == "NA"))
                .then(pl.col('Gene_2_id(3end_fusion_partner)'))
                .otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
            )
            
            # Base columns
            base_columns = [
                (gene1_expr + "::" + gene2_expr + '__' + pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("fusionTranscriptID"),
                (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("breakpointID"),
                (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2)).alias("5pStrand"),
                (pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2)).alias("3pStrand"),
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ]
            
            # FusionCatcher-specific columns
            fc_columns = []
            
            # Check if 'Predicted_effect' column exists and add it
            if 'Predicted_effect' in lazy_df.collect_schema().names():
                fc_columns.append(pl.col('Predicted_effect').alias('predictedEffect_FC'))
            else:
                fc_columns.append(pl.lit("NA").cast(pl.String).alias('predictedEffect_FC'))

            fc_columns.append(pl.col('Fusion_description').alias("fusionPairAnnotation_FC"))
            
            # Placeholder columns for other tools
            arr_columns = [
                pl.lit("NA").cast(pl.String).alias("5pSite_ARR"),
                pl.lit("NA").cast(pl.String).alias("3pSite_ARR"),
                pl.lit("NA").cast(pl.String).alias("mutationType_ARR"),
                pl.lit("NA").cast(pl.String).alias("confidenceLabel_ARR")
            ]
            
            sf_columns = [
                pl.lit("NA").cast(pl.String).alias("fusionPairAnnotation_SF"),
                pl.lit("NA").cast(pl.String).alias("junctionReadCount_SF"),
                pl.lit("NA").cast(pl.String).alias("spanningFragCount_SF"),
                pl.lit("NA").cast(pl.String).alias("largeAnchorSupport_SF")
            ]
            
            return lazy_df.select(base_columns + arr_columns + fc_columns + sf_columns)
        
        case 'STARFusion':
            # Extract gene names with fallback to gene IDs
            gene1_expr = (
                pl.when(pl.col('LeftGene').str.split("^").list.get(0) != "")
                .then(pl.col('LeftGene').str.split("^").list.get(0))
                .otherwise(
                    pl.col('LeftGene').str.split("^").list.get(1).str.split(".").list.get(0)
                )
            )
            gene2_expr = (
                pl.when(pl.col('RightGene').str.split("^").list.get(0) != "")
                .then(pl.col('RightGene').str.split("^").list.get(0))
                .otherwise(
                    pl.col('RightGene').str.split("^").list.get(1).str.split(".").list.get(0)
                )
            )
            
            # Handle breakpoints: Format from chr17:38243106:+ to 17:38243106
            left_breakpoint = (
                pl.col('LeftBreakpoint')
                .str.replace(r'^chr', '')  # Remove 'chr' prefix
                .str.split(':').list.slice(0, 2).list.join(':')  # Convert format
            )
            
            right_breakpoint = (
                pl.col('RightBreakpoint')
                .str.replace(r'^chr', '')  # Remove 'chr' prefix
                .str.split(':').list.slice(0, 2).list.join(':')  # Convert format
            )
            
            # Extract strands from breakpoint columns
            left_strand = pl.col('LeftBreakpoint').str.split(':').list.get(2)
            right_strand = pl.col('RightBreakpoint').str.split(':').list.get(2)
            
            # Base columns
            base_columns = [
                (gene1_expr + "::" + gene2_expr + '__' + left_breakpoint + "-" + right_breakpoint).alias("fusionTranscriptID"),
                (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                (left_breakpoint + "-" + right_breakpoint).alias("breakpointID"),
                left_strand.alias("5pStrand"),
                right_strand.alias("3pStrand"),
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ]
            
            # STARFusion-specific columns
            sf_columns = [
                pl.col('annots').cast(pl.String).alias("fusionPairAnnotation_SF"),
                pl.col('JunctionReadCount').cast(pl.String).alias("junctionReadCount_SF"),
                pl.col('SpanningFragCount').cast(pl.String).alias("spanningFragCount_SF"),
                pl.col('LargeAnchorSupport').cast(pl.String).alias("largeAnchorSupport_SF")
            ]
            
            # Placeholder columns for other tools
            arr_columns = [
                pl.lit("NA").cast(pl.String).alias("5pSite_ARR"),
                pl.lit("NA").cast(pl.String).alias("3pSite_ARR"),
                pl.lit("NA").cast(pl.String).alias("mutationType_ARR"),
                pl.lit("NA").cast(pl.String).alias("confidenceLabel_ARR")
            ]
            
            fc_columns = [
                pl.lit("NA").cast(pl.String).alias("predictedEffect_FC"),
                pl.lit("NA").cast(pl.String).alias("fusionPairAnnotation_FC")
            ]
            
            return lazy_df.select(base_columns + arr_columns + fc_columns + sf_columns)
            
        case _:
            raise ValueError(f"Unsupported tool name: {tool_name}")

def collate_fusion_data(
    sample_id: str, 
    output_filename: str, 
    input_files: List[Tuple[str, str]]
) -> None:
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

    for input_path, suffix in input_files:
        if not Path(input_path).exists():
            print(f"Error: File {input_path} does not exist. Skipping...")
            continue
        
        print('Setting tool name...')
        tool_name = tool_name_map.get(suffix)

        if not tool_name:
            print(f"Warning: Unknown tool suffix '{suffix}'. Skipping file {input_path}")
            continue

        # Extract sample number
        sample_num = extract_sample_num(input_path, suffix)
        if not sample_num:
            print(f"Error: Could not extract sample number from {input_path}.")
            sys.exit(1)
            
        print(f'Reading {tool_name} of {sample_id} TSV file...(sample ID: {sample_id})')

        # Create a lazy dataframe for each file
        try:
            lazy_df = wrangle_df(input_path, sample_id, sample_num, tool_name)
            # Append the lazy DataFrame to the list
            lazy_dfs.append(lazy_df)
        except Exception as e:
            print(f"Error processing {input_path}: {str(e)}")
            continue
    
    if not lazy_dfs:
        raise ValueError("No valid input files were processed. Aborting...")
    
    print(f"Concatenating lazy DataFrames from {len(lazy_dfs)} tools...")
        
    # Concatenate all lazy DataFrames
    combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)
    
    print("Concatenation completed. Collecting...")

    # Define categorical columns for the new structure
    base_categoricals = [
        "fusionTranscriptID",
        "fusionGenePair",
        "breakpointID",
        "5pStrand",
        "3pStrand",
        "originalTool",
        "sampleID"
    ]
    
    # Tool-specific categorical columns
    tool_categoricals = [
        # Arriba columns
        "5pSite_ARR",
        "3pSite_ARR",
        "mutationType_ARR",
        "confidenceLabel_ARR",
        # FusionCatcher columns
        "predictedEffect_FC",
        # STARFusion columns
        "largeAnchorSupport_SF"
    ]
    
    all_categoricals = base_categoricals + tool_categoricals
    
    # Integer columns
    ints = [
        "sampleNum"
    ]

    # Cast columns to appropriate types
    results = combined_lazy_df.with_columns(
        [pl.col(col).cast(pl.Categorical) for col in all_categoricals]
    ).with_columns(
        [pl.col(col).cast(pl.Int64) for col in ints]
    ).collect()
    
    print("DataFrame collected.")
    print(f"Final DataFrame shape: {results.shape}")
    print(f"Columns: {results.columns}")
    print(f"Saving as parquet and tsv files...")
        
    # Save as parquet and tsv
    results.write_parquet(f"{output_filename}.parquet")
    results.write_csv(f"{output_filename}.tsv", separator="\t")

    print("Done.")

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
                sys.exit(1)
                
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
            print(f"Error: {str(e)}")
            sys.exit(1)

if __name__ == "__main__":
    main()
    