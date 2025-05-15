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

# def natural_sort_key(s):
#     return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

def wrangle_df(file_path: str, sample_id:str, sample_num: str, tool_name: str) -> pl.LazyFrame:
    """Process input file based on the tool name and return a standardized lazy DataFrame."""
    lazy_df = pl.scan_csv(file_path, separator="\t")
    
    match tool_name:
        case 'Arriba':
            return lazy_df.select([
                (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
                (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
                (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
                (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("5pStrand"),
                (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("3pStrand"),
                (pl.col('site1')).alias("5pSite"),
				(pl.col('site2')).alias("3pSite"), 
                (pl.col('type')).alias("mutationType"),
				(pl.col('confidence')).alias("confidenceLabel"),
                # Add placeholder columns for STARFusion specific fields
                pl.lit('.').alias("largeAnchorSupport"),
                pl.lit(None).cast(pl.Int64).alias("junctionReadCount"),
                pl.lit(None).cast(pl.Int64).alias("spanningFragCount"),
                # Sample identification
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ])
        case 'FusionCatcher':
            # Handle NaN values in gene symbol columns by replacing with gene IDs
            gene1_expr = (
                pl.when(pl.col('Gene_1_symbol(5end_fusion_partner)').is_null() | (pl.col('Gene_1_symbol(5end_fusion_partner)') == ""))
                .then(pl.col('Gene_1_id(5end_fusion_partner)'))
                .otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
            )
            
            gene2_expr = (
                pl.when(pl.col('Gene_2_symbol(3end_fusion_partner)').is_null() | (pl.col('Gene_2_symbol(3end_fusion_partner)') == ""))
                .then(pl.col('Gene_2_id(3end_fusion_partner)'))
                .otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
            )
            
            base_columns = [
                (gene1_expr + "::" + gene2_expr + '__' + pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("fusionTranscriptID"),
                (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':') + "-" + pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(':').list.slice(0, 2).list.join(':')).alias("breakpointID"),
                (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2)).alias("5pStrand"),
                (pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2)).alias("3pStrand"),
            ]
            
            # Check if 'Predicted_effect' column exists
            if 'Predicted_effect' in lazy_df.collect_schema().names():
                predicted_effect_columns = [
                    pl.col('Predicted_effect').str.extract(r'^([^/]+)(?:/|$)').alias('5pSite'),
                    pl.when(pl.col('Predicted_effect').str.contains('/'))
                        .then(pl.col('Predicted_effect').str.extract(r'/(.+)$'))
                        .otherwise(pl.lit('.'))
                        .alias('3pSite')
                    ]
            else:
                predicted_effect_columns = [
                    pl.lit('.').alias("5pSite"),
                    pl.lit('.').alias("3pSite")
                    ]

            additional_columns = [
                pl.lit('.').alias("mutationType"),
                pl.lit('.').alias("confidenceLabel"),
                # Add placeholder columns for STARFusion specific fields
                pl.lit('.').alias("largeAnchorSupport"),
                pl.lit(None).cast(pl.Int64).alias("junctionReadCount"),
                pl.lit(None).cast(pl.Int64).alias("spanningFragCount"),
                # Sample identification
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ]
            
            return lazy_df.select(base_columns + predicted_effect_columns + additional_columns)
        
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
            
            return lazy_df.select([
                # Core fusion identification columns
                (gene1_expr + "::" + gene2_expr + '__' + left_breakpoint + "-" + right_breakpoint).alias("fusionTranscriptID"),
                (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                (left_breakpoint + "-" + right_breakpoint).alias("breakpointID"),
                left_strand.alias("5pStrand"),
                right_strand.alias("3pStrand"),
                pl.lit('.').alias("5pSite"),
                pl.lit('.').alias("3pSite"),
                pl.lit('.').alias("mutationType"),
                pl.lit('.').alias("confidenceLabel"),
                # STARFusion specific columns
                pl.col('LargeAnchorSupport').alias("largeAnchorSupport"),
                pl.col('JunctionReadCount').alias("junctionReadCount"),
                pl.col('SpanningFragCount').alias("spanningFragCount"),
                # Sample identification
                pl.lit(tool_name).alias("originalTool"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_num).cast(pl.Int64).alias("sampleNum"),
                pl.lit(sample_num).cast(pl.Utf8).str.zfill(4).alias("sampleNum_Padded")
            ])
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

    # Define categorical columns including the new STARFusion specific columns
    categoricals = [
        "fusionTranscriptID",
        "fusionGenePair",
        "breakpointID",
        "5pStrand",
        "3pStrand",
        "5pSite",
        "3pSite",
        "mutationType",
        "confidenceLabel",
        # Additional STARFusion specific columns
        "largeAnchorSupport",
        # Sample identification columns
        "originalTool",
        "sampleID"
    ]

    ints = [
        "junctionReadCount",
        "spanningFragCount",
        "sampleNum"
    ]

    # Cast columns to appropriate types
    results = combined_lazy_df.with_columns(
        [pl.col(col).cast(pl.Categorical) for col in categoricals]
    ).with_columns(
        [pl.col(col).cast(pl.Int64) for col in ints]
    ).collect()
    
    print("DataFrame collected.")
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
            input_files.append((os.path.abspath(args.arriba), 'arr'))
        if args.fusioncatcher:
            input_files.append((os.path.abspath(args.fusioncatcher), 'fc'))
        if args.starfusion:
            input_files.append((os.path.abspath(args.starfusion), 'sf'))
            
        # Add additional inputs if provided
        if args.inputs:
            if len(args.inputs) % 2 != 0:
                print("Error: Additional inputs must be provided as pairs of file path and tool suffix.")
                sys.exit(1)
                
            for i in range(0, len(args.inputs), 2):
                input_files.append((os.path.abspath(args.inputs[i]), args.inputs[i+1]))
        
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
    
    else:
        # Legacy positional argument handling for backward compatibility
        if len(sys.argv) < 4:
            print("Usage: collate-FTs--nf.py <sample_id> <output_file> [<file_path> <tool_suffix>]...")
            print("Example: collate-FTs--nf.py 123T output path/to/123T_arr.tsv arr path/to/123T_fc.tsv fc path/to/123T_sf.tsv sf")
            print("\nAlternatively, use the new argument format:")
            print("collate-FTs--nf.py --sample_id 123T --output output --arriba path/to/123T_arr.tsv --fusioncatcher path/to/123T_fc.tsv")
            sys.exit(1)
            
        sample_id = sys.argv[1]
        output_file = sys.argv[2]
        
        # Check if the number of remaining arguments is even (each file needs a suffix)
        if (len(sys.argv) - 3) % 2 != 0:
            print("Error: Each input file provided must be tagged with its corresponding tool suffix as an argument.")
            sys.exit(1)
        
        # Create a list of tuples (file_path, suffix)
        input_files = []
        for i in range(3, len(sys.argv), 2):
            if i+1 < len(sys.argv):
                input_files.append((os.path.abspath(sys.argv[i]), sys.argv[i+1]))
        
        # Call the collate function
        try:
            collate_fusion_data(
                sample_id=sample_id,
                output_filename=output_file,
                input_files=input_files
            )
        except Exception as e:
            print(f"Error: {str(e)}")
            sys.exit(1)


if __name__ == "__main__":
    main()
