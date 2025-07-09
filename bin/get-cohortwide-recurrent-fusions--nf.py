#!/usr/bin/env python3
"""
Script to filter for recurrent fusions based on frequency threshold
"""
import sys
import math
import argparse
from pathlib import Path
import polars as pl

def filter_recurrent_fusions(input_file, threshold, output_file, report_file):
    """
    Filter fusions based on recurrence frequency across samples using breakpointID
    """
    try:
        # Read the cohort TSV file
        print(f"Reading cohort file: {input_file}")
        df = pl.read_csv(input_file, separator='\t')
        print(f"Total rows in cohort file: {len(df)}")
        
        # Get total number of samples directly from the data
        total_samples = df['sampleID'].n_unique()
        print(f"Total samples in cohort: {total_samples}")
        
        # Calculate minimum samples needed based on threshold
        min_samples_required = int(math.ceil(threshold * total_samples)) 
        # if min_samples_required == 1, we have to override the threshold so that min_samples_required is at least 2
        if min_samples_required <= 1:
            min_samples_required = 2
            print("Warning: Total samples in the input cohort is too low for the specified threshold. Setting minimum number of samples to determine recurrence to 2.")
            # compute and update the threshold accordingly
            threshold = min_samples_required / total_samples
            print(f"Updated threshold to {threshold*100:.1f}% to ensure at least 2 samples are considered for fusion recurrence.")
        else:
            print(f"Minimum samples for recurrence at {threshold*100:.1f}% threshold: {min_samples_required}")
        
        # Calculate fusion frequencies using breakpointID for granularity
        print("Calculating fusion frequencies by breakpointID...")
        
        # Group by breakpointID and count unique samples
        fusion_freq = (df
                      .group_by('breakpointID')
                      .agg([
                          pl.col('sampleID').n_unique().alias('sample_count'),
                          pl.col('sampleID').unique().alias('samples_with_fusion'),
                          pl.col('fusionGenePair').first().alias('fusionGenePair'),  # Keep gene pair for reference
                          pl.len().alias('total_occurrences')
                      ]))
        
        # Calculate frequency percentage
        fusion_freq = fusion_freq.with_columns([
            (pl.col('sample_count') / total_samples * 100).alias('frequency_percent')
        ])
        
        # Sort by frequency
        fusion_freq = fusion_freq.sort('frequency_percent', descending=True)
        
        print(f"Found {len(fusion_freq)} unique breakpoints")
        
        # Apply filters
        threshold_percent = threshold * 100
        print(f"Applying filters:")
        print(f"  - Frequency for recurrence: {threshold_percent:.2f}%")
        print(f"  - Minimum samples: {min_samples_required}")
        
        recurrent_fusions = fusion_freq.filter(
            (pl.col('frequency_percent') >= threshold_percent) &
            (pl.col('sample_count') >= min_samples_required)
        )
        
        print(f"Recurrent breakpoints found: {len(recurrent_fusions)}")
        
        if len(recurrent_fusions) == 0:
            print("Warning: No recurrent fusions found with current criteria!")
            # Create empty output file with headers
            empty_df = df.head(0)
            empty_df.write_csv(output_file, separator='\t')
        else:
            # Filter original data to keep only recurrent fusions
            recurrent_breakpoints = recurrent_fusions['breakpointID'].to_list()
            
            filtered_df = df.filter(pl.col('breakpointID').is_in(recurrent_breakpoints))
            
            # Add frequency rank for sorting
            freq_rank = (recurrent_fusions
                        .with_row_count('frequency_rank', offset=1)
                        .select(['breakpointID', 'frequency_rank']))
            
            filtered_df = filtered_df.join(freq_rank, on='breakpointID', how='left')
            filtered_df = filtered_df.sort(['frequency_rank', 'sampleID'])
            filtered_df = filtered_df.drop('frequency_rank')
            
            # Write filtered results
            filtered_df.write_csv(output_file, separator='\t')
            print(f"Filtered data written to: {output_file}")
            print(f"Rows in filtered file: {len(filtered_df)}")
        
        # Write frequency report
        with open(report_file, 'w') as f:
            f.write("Fusion Frequency Report\n")
            f.write("=====================\n\n")
            f.write(f"Total samples: {total_samples}\n")
            f.write(f"Total unique breakpoints: {len(fusion_freq)}\n")
            f.write(f"Recurrence threshold: {threshold_percent:.2f}%\n")
            f.write(f"Minimum samples required: {min_samples_required}\n")
            f.write(f"Recurrent breakpoints found: {len(recurrent_fusions)}\n\n")
            
            if len(recurrent_fusions) > 0:
                f.write(f"\nRecurrent Breakpoints (>= {threshold_percent:.2f}%):\n")
                f.write("-" * 100 + "\n")
                # Write header with fixed-width columns
                f.write(f"{'BreakpointID':<40} {'FusionGenePair':<40} {'SampleCount':<6} {'Frequency%':<8} {'TotalOccurrences':<6} {'Samples'}\n")
                f.write("-" * 130 + "\n")
                for row in recurrent_fusions.iter_rows():
                    breakpoint_id, sample_count, samples, gene_pair, total_occ, freq_pct = row
                    sample_list = ",".join(samples[:5])  # Show first 5 samples
                    if len(samples) > 5:
                        sample_list += f"... (+{len(samples)-5} more)"
                    f.write(f"{breakpoint_id:<35} {gene_pair:<40} {sample_count:<6} {freq_pct:<8.2f} {total_occ:<6} {sample_list}\n")
        
        print(f"Frequency report written to: {report_file}")
        
    except Exception as e:
        print(f"Error filtering fusions: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Filter for recurrent fusions')
    parser.add_argument('--input', required=True,
                       help='Input cohort TSV file')
    parser.add_argument('--threshold', type=float, default=0.005,
                       help='Frequency threshold (default: 0.005 = 0.5%%)')
    parser.add_argument('--output', required=True,
                       help='Output filtered TSV file')
    parser.add_argument('--report', required=True,
                       help='Output frequency report file')
    
    args = parser.parse_args()
    
    # Validate threshold
    if args.threshold <= 0 or args.threshold > 1:
        print("Error: Threshold must be between 0 and 1 (e.g., 0.005 for 0.5%)")
        sys.exit(1)
    
    # Filter for recurrent fusions
    filter_recurrent_fusions(
        args.input, 
        args.threshold, 
        args.output, 
        args.report
    )
    
    print(f"\nFiltering completed successfully!")

if __name__ == "__main__":
    main()
    
