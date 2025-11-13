#!/usr/bin/env python3
"""
Script to visualize recurrent fusion distributions
"""
import sys
import argparse
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

def create_fusion_visualization(unthresholded_report, output_html, min_cohort_size=10, report_type="Unthresholded"):
    """
    Create interactive HTML visualizations of recurrent fusion distributions
    
    Args:
        unthresholded_report: Path to fusion report TSV file
        output_html: Output HTML file path
        min_cohort_size: Minimum cohort size to generate visualization
        report_type: Type of report ("Unthresholded" or "Thresholded")
    """
    try:
        # Read the fusion report TSV
        print(f"Reading fusion report: {unthresholded_report}")
        df = pl.read_csv(unthresholded_report, separator='\t')
        
        if len(df) == 0:
            print("Warning: No recurrent fusions found in report. Skipping visualization.")
            return False
        
        print(f"Found {len(df)} recurrent fusion breakpoints")
        
        # Get cohort size from the report (count unique samples across all fusions)
        all_samples = set()
        for samples_str in df['SamplesWithFusion']:
            all_samples.update(samples_str.split(','))
        cohort_size = len(all_samples)
        
        print(f"Cohort size (fusion-positive samples): {cohort_size}")
        
        if cohort_size < min_cohort_size:
            print(f"Cohort size ({cohort_size}) is below minimum threshold ({min_cohort_size}). Skipping visualization.")
            return False
        
        # Create subplots: pie chart and bar chart
        subtitle_prefix = "All " if report_type == "Unthresholded" else "Thresholded "
        min_samples_text = "(≥2 samples)" if report_type == "Unthresholded" else ""
        
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(
                f'Distribution of {subtitle_prefix}Recurrent Fusion Gene Pairs',
                f'Top 20 {subtitle_prefix}Recurrent Fusions by Sample Count'
            ),
            specs=[[{"type": "pie"}], [{"type": "bar"}]],
            vertical_spacing=0.15,
            row_heights=[0.5, 0.5]
        )
        
        # Aggregate by gene pair for pie chart
        gene_pair_counts = (df
                           .group_by('fusionGenePair')
                           .agg([
                               pl.col('SampleCount').sum().alias('total_samples'),
                               pl.len().alias('num_breakpoints')
                           ])
                           .sort('total_samples', descending=True))
        
        # For pie chart: show top 10, group rest as "Other"
        top_n = 10
        if len(gene_pair_counts) > top_n:
            top_pairs = gene_pair_counts.head(top_n)
            other_count = gene_pair_counts.tail(len(gene_pair_counts) - top_n)['total_samples'].sum()
            
            labels = top_pairs['fusionGenePair'].to_list() + ['Other']
            values = top_pairs['total_samples'].to_list() + [other_count]
        else:
            labels = gene_pair_counts['fusionGenePair'].to_list()
            values = gene_pair_counts['total_samples'].to_list()
        
        # Add pie chart
        fig.add_trace(
            go.Pie(
                labels=labels,
                values=values,
                textinfo='label+percent',
                hovertemplate='<b>%{label}</b><br>' +
                             'Total Samples: %{value}<br>' +
                             'Percentage: %{percent}<br>' +
                             '<extra></extra>',
                marker=dict(line=dict(color='white', width=2))
            ),
            row=1, col=1
        )
        
        # Prepare data for bar chart (top 20 individual breakpoints)
        top_breakpoints = df.sort('SampleCount', descending=True).head(20)
        
        # Create hover text with sample IDs (show first 5)
        hover_text = []
        for row in top_breakpoints.iter_rows(named=True):
            samples = row['SamplesWithFusion'].split(',')[:5]
            sample_preview = ', '.join(samples)
            if row['SampleCount'] > 5:
                sample_preview += f" ... (+{row['SampleCount'] - 5} more)"
            hover_text.append(
                f"<b>{row['fusionGenePair']}</b><br>" +
                f"Breakpoint: {row['breakpointID']}<br>" +
                f"Samples: {row['SampleCount']} ({row['FrequencyPercent']:.2f}%)<br>" +
                f"Occurrences: {row['TotalOccurrences']}<br>" +
                f"Found in: {sample_preview}"
            )
        
        # Add bar chart
        fig.add_trace(
            go.Bar(
                x=top_breakpoints['SampleCount'].to_list(),
                y=[f"{row['fusionGenePair']}<br>({row['breakpointID'][:20]}...)" 
                   for row in top_breakpoints.iter_rows(named=True)],
                orientation='h',
                marker=dict(
                    color=top_breakpoints['FrequencyPercent'].to_list(),
                    colorscale='Viridis',
                    colorbar=dict(
                        title="Frequency %",
                        x=1.15
                    )
                ),
                hovertext=hover_text,
                hoverinfo='text',
                text=top_breakpoints['SampleCount'].to_list(),
                textposition='outside'
            ),
            row=2, col=1
        )
        
        # Update layout
        title_text = f'{report_type} Recurrent Fusion Distribution Analysis<br>'
        if report_type == "Unthresholded":
            subtitle_text = f'<sub>Cohort: {cohort_size} fusion-positive samples | ' + \
                          f'{len(df)} recurrent breakpoints (≥2 samples)</sub>'
        else:
            subtitle_text = f'<sub>Cohort: {cohort_size} fusion-positive samples | ' + \
                          f'{len(df)} significant recurrent breakpoints</sub>'
        
        fig.update_layout(
            title=dict(
                text=title_text + subtitle_text,
                x=0.5,
                xanchor='center'
            ),
            showlegend=False,
            height=1000,
            template='plotly_white',
            font=dict(size=12)
        )
        
        # Update axes for bar chart
        fig.update_xaxes(title_text="Number of Samples", row=2, col=1)
        fig.update_yaxes(title_text="Fusion (Gene Pair & Breakpoint)", row=2, col=1)
        
        # Save to HTML
        fig.write_html(
            output_html,
            include_plotlyjs='cdn',
            config={
                'displayModeBar': True,
                'displaylogo': False,
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': 'recurrent_fusions_visualization',
                    'height': 1000,
                    'width': 1200,
                    'scale': 2
                }
            }
        )
        
        print(f"Visualization written to: {output_html}")
        return True
        
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description='Visualize recurrent fusion distributions')
    parser.add_argument('--input', required=True,
                       help='Input fusion report TSV file')
    parser.add_argument('--output', required=True,
                       help='Output HTML file for visualization')
    parser.add_argument('--min-cohort-size', type=int, default=10,
                       help='Minimum cohort size to generate visualization (default: 10)')
    parser.add_argument('--report-type', choices=['Unthresholded', 'Thresholded'], 
                       default='Unthresholded',
                       help='Type of report being visualized (default: Unthresholded)')
    
    args = parser.parse_args()
    
    # Validate min_cohort_size
    if args.min_cohort_size < 1:
        print("Error: --min-cohort-size must be at least 1")
        sys.exit(1)
    
    # Create visualization
    success = create_fusion_visualization(
        args.input,
        args.output,
        args.min_cohort_size,
        args.report_type
    )
    
    if success:
        print("\nVisualization completed successfully!")
    else:
        print("\nVisualization skipped or failed.")
        sys.exit(0)  # Don't fail the pipeline, just skip

if __name__ == "__main__":
    main()
	