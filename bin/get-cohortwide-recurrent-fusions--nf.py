#!/usr/bin/env python3
"""
Script to filter for recurrent fusions based on frequency threshold
"""
import sys
import math
import argparse
from pathlib import Path
import polars as pl

def filter_recurrent_fusions(input_file, threshold, output_file, report_file, total_samples, unthresholded_output=None):
    """
    Filter fusions based on recurrence frequency across samples using breakpointID
    """
    try:
        # Read the cohort TSV file
        print(f"Reading cohort file: {input_file}")
        df = pl.read_csv(input_file, separator='\t')
        print(f"Total rows in cohort file: {len(df)}")
        
        # Get samples with fusions from file
        samples_with_fusions = df['sampleID'].n_unique()
        print(f"Total input samples [cohort count]: {total_samples}")
        print(f"Samples with (validated) fusions [fusion-positives]: {samples_with_fusions}")
        
        # Calculate minimum samples needed based on threshold
        # Use fusion-positive samples as denominator for biologically meaningful recurrence
        min_samples_required = int(math.ceil(threshold * samples_with_fusions)) 
        # if min_samples_required == 1, we have to override the threshold so that min_samples_required is at least 2
        if min_samples_required <= 1:
            min_samples_required = 2
            print("Warning: Total fusion-positive samples in the input cohort is too low for the specified threshold. Setting minimum number of samples to determine recurrence to 2.")
            # compute and update the threshold accordingly
            threshold = min_samples_required / samples_with_fusions
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
        
        # Calculate frequency percentage based on fusion-positive samples
        fusion_freq = fusion_freq.with_columns([
            (pl.col('sample_count') / samples_with_fusions * 100).alias('frequency_percent')
        ])
        
        # Sort by frequency
        fusion_freq = fusion_freq.sort('frequency_percent', descending=True)
        
        print(f"Found {len(fusion_freq)} unique breakpoints")
        
        # Filter for unthresholded recurrent fusions (appearing in >= 2 samples)
        if unthresholded_output:
            print("\nGenerating unthresholded recurrent fusions (>= 2 samples)...")
            unthresholded_recurrent = fusion_freq.filter(pl.col('sample_count') >= 2)
            print(f"Unthresholded recurrent breakpoints (>= 2 samples): {len(unthresholded_recurrent)}")
            
            if len(unthresholded_recurrent) > 0:
                # Filter original data for unthresholded recurrent fusions
                unthresholded_breakpoints = unthresholded_recurrent['breakpointID'].to_list()
                unthresholded_df = df.filter(pl.col('breakpointID').is_in(unthresholded_breakpoints))
                
                # Add frequency rank for sorting
                unthresh_freq_rank = (unthresholded_recurrent
                                     .with_row_index('frequency_rank', offset=1)
                                     .select(['breakpointID', 'frequency_rank']))
                
                unthresholded_df = unthresholded_df.join(unthresh_freq_rank, on='breakpointID', how='left')
                unthresholded_df = unthresholded_df.sort(['frequency_rank', 'sampleID'])
                unthresholded_df = unthresholded_df.drop('frequency_rank')
                
                # Write unthresholded results
                unthresholded_df.write_csv(unthresholded_output, separator='\t')
                print(f"Unthresholded recurrent fusions written to: {unthresholded_output}")
                print(f"Rows in unthresholded file: {len(unthresholded_df)}")
                
                # Write unthresholded human-readable report
                unthresh_txt_report = str(Path(unthresholded_output).parent / (Path(unthresholded_output).stem + '.freq-report.txt'))
                with open(unthresh_txt_report, 'w') as f:
                    f.write("Unthresholded Recurrent Fusion Frequency Report\n")
                    f.write("=" * 100 + "\n\n")
                    f.write(f"Total input samples [cohort count]: {total_samples}\n")
                    f.write(f"Samples with (validated) fusions [fusion-positives]: {samples_with_fusions}\n")
                    f.write(f"Samples without fusions: {total_samples - samples_with_fusions}\n\n")
                    f.write(f"NOTE: This report shows ALL recurrent fusions (appearing in >= 2 samples)\n")
                    f.write(f"      without applying frequency thresholding.\n\n")
                    f.write(f"Total unique breakpoints: {len(fusion_freq)}\n")
                    f.write(f"Unthresholded recurrent breakpoints (>= 2 samples): {len(unthresholded_recurrent)}\n\n")
                    f.write("-" * 100 + "\n")
                    f.write(f"{'BreakpointID':<40} {'FusionGenePair':<40} {'SampleCount':<6} {'Frequency%':<8} {'TotalOccurrences':<6} {'Samples'}\n")
                    f.write("-" * 130 + "\n")
                    for row in unthresholded_recurrent.iter_rows():
                        breakpoint_id, sample_count, samples, gene_pair, total_occ, freq_pct = row
                        sample_list = ",".join(samples[:5])  # Show first 5 samples
                        if len(samples) > 5:
                            sample_list += f"... (+{len(samples)-5} more)"
                        f.write(f"{breakpoint_id:<35} {gene_pair:<40} {sample_count:<6} {freq_pct:<8.2f} {total_occ:<6} {sample_list}\n")
                
                print(f"Unthresholded human-readable report written to: {unthresh_txt_report}")
                
                # Write unthresholded TSV report
                unthresh_report = str(Path(unthresholded_output).with_suffix('.freq-report.tsv'))
                unthresh_tsv_data = unthresholded_recurrent.with_columns([
                    pl.col('samples_with_fusion').list.join(',').alias('samples_list')
                ]).select([
                    'fusionGenePair',
                    'breakpointID', 
                    'sample_count',
                    'frequency_percent',
                    'total_occurrences',
                    'samples_list'
                ]).rename({
                    'sample_count': 'SampleCount',
                    'frequency_percent': 'FrequencyPercent',
                    'total_occurrences': 'TotalOccurrences',
                    'samples_list': 'SamplesWithFusion'
                }).sort(['fusionGenePair', 'breakpointID'])
                
                unthresh_tsv_data.write_csv(unthresh_report, separator='\t')
                print(f"Unthresholded frequency report written to: {unthresh_report}")
            else:
                print("Warning: No breakpoints found in >= 2 samples!")
                empty_df = df.head(0)
                empty_df.write_csv(unthresholded_output, separator='\t')

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
            # Create empty output file
            empty_df = df.head(0)
            empty_df.write_csv(output_file, separator='\t')
        else:
            # Filter original data to keep only recurrent fusions
            recurrent_breakpoints = recurrent_fusions['breakpointID'].to_list()
            
            filtered_df = df.filter(pl.col('breakpointID').is_in(recurrent_breakpoints))
            
            # Add frequency rank for sorting
            freq_rank = (recurrent_fusions
                        .with_row_index('frequency_rank', offset=1)
                        .select(['breakpointID', 'frequency_rank']))
            
            filtered_df = filtered_df.join(freq_rank, on='breakpointID', how='left')
            filtered_df = filtered_df.sort(['frequency_rank', 'sampleID'])
            filtered_df = filtered_df.drop('frequency_rank')
            
            # Write filtered results
            filtered_df.write_csv(output_file, separator='\t')
            print(f"Filtered data written to: {output_file}")
            print(f"Rows in filtered file: {len(filtered_df)}")
        
        # Write frequency report (human-readable)
        with open(report_file, 'w') as f:
            f.write("Fusion Frequency Report\n")
            f.write("=====================\n\n")
            f.write(f"Total input samples [cohort count]: {total_samples}\n")
            f.write(f"Samples with (validated) fusions [fusion-positives]: {samples_with_fusions}\n")
            f.write(f"Samples without fusions: {total_samples - samples_with_fusions}\n\n")
            f.write(f"NOTE: Recurrence calculated as proportion of fusion-positive samples ({samples_with_fusions})\n")
            f.write(f"      to ensure biologically meaningful frequency estimates.\n\n")
            f.write(f"Total unique breakpoints: {len(fusion_freq)}\n")
            
            if unthresholded_output:
                unthresholded_count = len(fusion_freq.filter(pl.col('sample_count') >= 2))
                f.write(f"Unthresholded recurrent breakpoints (>= 2 samples): {unthresholded_count}\n")
            
            f.write(f"Recurrence threshold: {threshold_percent:.2f}% of fusion-positive samples\n")
            f.write(f"Minimum samples required: {min_samples_required}\n")
            f.write(f"Recurrent breakpoints retained after applying recurrence threshold: {len(recurrent_fusions)}\n\n")
            
            if len(recurrent_fusions) > 0:
                f.write(f"\nRecurrent Breakpoints (>= {threshold_percent:.2f}% of fusion-positive samples):\n")
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
        
        # Write frequency report as TSV (machine-readable)
        tsv_report_file = str(Path(report_file).parent / (Path(report_file).stem + '.tsv'))
        
        # Prepare TSV data with samples as comma-separated string
        tsv_data = recurrent_fusions.with_columns([
            pl.col('samples_with_fusion').list.join(',').alias('samples_list')
        ]).select([
            'fusionGenePair',
            'breakpointID', 
            'sample_count',
            'frequency_percent',
            'total_occurrences',
            'samples_list'
        ]).rename({
            'sample_count': 'SampleCount',
            'frequency_percent': 'FrequencyPercent',
            'total_occurrences': 'TotalOccurrences',
            'samples_list': 'SamplesWithFusion'
        }).sort(['fusionGenePair', 'breakpointID'])
        
        if len(recurrent_fusions) > 0:
            tsv_data.write_csv(tsv_report_file, separator='\t')
            print(f"TSV frequency report written to: {tsv_report_file}")
        else:
            # Create empty TSV with headers
            empty_tsv = pl.DataFrame({
                'fusionGenePair': [],
                'breakpointID': [],
                'SampleCount': [],
                'FrequencyPercent': [],
                'TotalOccurrences': [],
                'SamplesWithFusion': []
            }, schema={
                'fusionGenePair': pl.Utf8,
                'breakpointID': pl.Utf8,
                'SampleCount': pl.Int64,
                'FrequencyPercent': pl.Float64,
                'TotalOccurrences': pl.Int64,
                'SamplesWithFusion': pl.Utf8
            })
            empty_tsv.write_csv(tsv_report_file, separator='\t')
            print(f"Empty TSV frequency report written to: {tsv_report_file}")
        
        print(f"Human-readable frequency report written to: {report_file}")
        
    except Exception as e:
        print(f"Error filtering fusions: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Filter for recurrent fusions')
    parser.add_argument('--input', required=True,
                       help='Input cohort TSV file')
    parser.add_argument('--total-samples', type=int, required=True,
                       help='Total number of samples in cohort (including fusion-negative samples)')
    parser.add_argument('--threshold', type=float, default=0.005,
                       help='Frequency threshold (default: 0.005 = 0.5%% of fusion-positive samples)')
    parser.add_argument('--output', required=True,
                       help='Output filtered TSV file')
    parser.add_argument('--report', required=True,
                       help='Output frequency report file')
    parser.add_argument('--unthresholded-output', required=False, default='unthresholded_recurrent_fusions.tsv',
                       help='Output file for unthresholded recurrent fusions (>= 2 samples) (default: unthresholded_recurrent_fusions.tsv). Set to empty string to disable.')
    
    args = parser.parse_args()
    
    # Validate total_samples
    if args.total_samples <= 0:
        print("Error: --total-samples must be a positive integer")
        sys.exit(1)
    
    # Validate threshold
    if args.threshold <= 0 or args.threshold > 1:
        print("Error: Threshold must be between 0 and 1 (e.g., 0.005 for 0.5%)")
        sys.exit(1)
    
    # Filter for recurrent fusions
    filter_recurrent_fusions(
        args.input, 
        args.threshold, 
        args.output, 
        args.report,
        args.total_samples,
        args.unthresholded_output if args.unthresholded_output else None
    )
    
    print(f"\nFiltering completed successfully!")

if __name__ == "__main__":
    main()
