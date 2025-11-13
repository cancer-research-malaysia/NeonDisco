#!/usr/bin/env python3
"""
Script to create interactive heatmap of patient-fusion coverage
"""
import sys
import argparse
import polars as pl
import plotly.graph_objects as go
import plotly.figure_factory as ff
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

def create_fusion_heatmap(report_tsv, output_html, max_fusions=None, report_type="Unthresholded", cluster=True):
    """
    Create interactive heatmap showing which patients have which fusions
    
    Args:
        report_tsv: Path to fusion report TSV file
        output_html: Output HTML file path
        max_fusions: Maximum number of fusions to display (None = all)
        report_type: Type of report ("Unthresholded" or "Thresholded")
        cluster: Whether to perform hierarchical clustering (default: True)
    """
    try:
        # Read the fusion report TSV
        print(f"Reading fusion report: {report_tsv}")
        df = pl.read_csv(report_tsv, separator='\t')
        
        if len(df) == 0:
            print("Warning: No recurrent fusions found in report. Skipping heatmap.")
            return False
        
        print(f"Found {len(df)} recurrent fusion breakpoints")
        
        # Limit number of fusions if specified
        if max_fusions and len(df) > max_fusions:
            print(f"Limiting to top {max_fusions} fusions by sample count")
            df = df.sort('SampleCount', descending=True).head(max_fusions)
        
        # Create fusion labels (gene pair + truncated breakpoint)
        fusion_labels = []
        for row in df.iter_rows(named=True):
            # Check if row is available before accessing keys (just in case polars is empty)
            if row:
                bp_short = row['breakpointID'][:30] + "..." if len(row['breakpointID']) > 30 else row['breakpointID']
                label = f"{row['fusionGenePair']}<br><sub>{bp_short}</sub>"
                fusion_labels.append(label)
        
        # Get all unique patients across all fusions
        all_patients = set()
        for samples_str in df['SamplesWithFusion']:
            all_patients.update(samples_str.split(','))
        all_patients = sorted(list(all_patients))
        
        print(f"Total patients with fusions: {len(all_patients)}")
        
        # Create binary matrix: rows = fusions, columns = patients
        matrix = np.zeros((len(df), len(all_patients)), dtype=int)
        patient_to_idx = {patient: idx for idx, patient in enumerate(all_patients)}
        
        for fusion_idx, row in enumerate(df.iter_rows(named=True)):
            patients = row['SamplesWithFusion'].split(',')
            for patient in patients:
                # Check if patient exists in the map (should always be true)
                if patient in patient_to_idx:
                    patient_idx = patient_to_idx[patient]
                    matrix[fusion_idx, patient_idx] = 1
        
        # Calculate patient fusion burden (how many fusions each patient has)
        patient_fusion_counts = matrix.sum(axis=0)

        # --- Clustering / Sorting Logic (Fixed) ---
        
        # Initialize default order variables
        fusion_order = list(range(len(df)))
        patient_order = list(range(len(all_patients)))
        
        # Perform hierarchical clustering if requested and we have enough data
        if cluster and len(df) >= 2 and len(all_patients) >= 2:
            print("Performing hierarchical clustering...")
            
            # Cluster rows (fusions) - using Jaccard distance for binary data
            # Only cluster fusions if there are at least two
            fusion_distances = pdist(matrix, metric='jaccard')
            fusion_linkage = linkage(fusion_distances, method='average')
            fusion_dendrogram = dendrogram(fusion_linkage, no_plot=True)
            fusion_order = fusion_dendrogram['leaves']
            
            # Cluster columns (patients) - using Jaccard distance for binary data
            # Only cluster patients if there are at least two
            patient_distances = pdist(matrix.T, metric='jaccard')
            patient_linkage = linkage(patient_distances, method='average')
            patient_dendrogram = dendrogram(patient_linkage, no_plot=True)
            patient_order = patient_dendrogram['leaves']

            print(f"Clustering complete: {len(df)} fusions × {len(all_patients)} patients")

        else:
            # Sort patients by fusion burden (descending) if not clustering
            print("Sorting by fusion burden (not clustering)...")
            # No fusion reordering (fusion_order remains 0, 1, 2...)
            patient_order = np.argsort(-patient_fusion_counts)
        
        # Apply final order to all necessary data structures
        matrix_final = matrix[fusion_order, :][:, patient_order]
        fusion_labels_final = [fusion_labels[i] for i in fusion_order]
        patients_final = [all_patients[i] for i in patient_order]
        patient_fusion_counts_final = patient_fusion_counts[patient_order]
        
        # Create a new Polars DataFrame reordered by fusion_order for hover text alignment
        df_final = df[fusion_order]
        # --- Heatmap and Hover Text Creation (Fixed) ---
        # Create hover text with detailed info
        hover_text = []
        # Iterate over the final, ordered list of fusions (now in df_final)
        for row in df_final.iter_rows(named=True):
            hover_row = []
            patients_with_fusion = set(row['SamplesWithFusion'].split(','))
            # Iterate over the final, ordered list of patients
            for patient in patients_final:
                if patient in patients_with_fusion:
                    text = (f"<b>{row['fusionGenePair']}</b><br>"
                           f"Patient: {patient}<br>"
                           f"Breakpoint: {row['breakpointID']}<br>"
                           f"Frequency: {row['FrequencyPercent']:.2f}%<br>"
                           f"Total samples with fusion: {row['SampleCount']}")
                else:
                    text = f"Patient: {patient}<br>No fusion detected"
                hover_row.append(text)
            hover_text.append(hover_row)
        
        # Create the heatmap
        fig = go.Figure(data=go.Heatmap(
            z=matrix_final, # Use the final, sorted/clustered matrix
            x=patients_final, # Use the final, sorted/clustered patient list
            y=fusion_labels_final, # Use the final, sorted/clustered fusion labels
            hovertext=hover_text,
            hoverinfo='text',
            colorscale=[[0, '#f0f0f0'], [1, '#d62728']],  # Gray for absent, red for present
            showscale=False,
            xgap=1,
            ygap=1
        ))
        
        # Add annotation showing patient fusion burden
        patient_burden_text = [f"{int(count)} fusions" if count > 0 else "0" 
                              for count in patient_fusion_counts_final]
        
        # Determine appropriate height based on number of fusions
        base_height = 400
        height_per_fusion = 25
        total_height = max(base_height, base_height + len(df) * height_per_fusion)
        
        # Update layout
        title_text = f'{report_type} Recurrent Fusion Coverage Matrix'
        subtitle_text = (f'<sub>{len(all_patients)} patients × {len(df)} fusions | '
                        f'Red = fusion present, Gray = fusion absent</sub>')
        
        fig.update_layout(
            title=dict(
                text=f'{title_text}<br>{subtitle_text}',
                x=0.5,
                xanchor='center',
                font=dict(size=16)
            ),
            xaxis=dict(
                title='Patient ID',
                tickangle=-90,
                tickfont=dict(size=8),
                side='bottom',
                # Add patient fusion burden as a secondary x-axis label (optional, can be complex in Plotly)
                # For simplicity, we stick to the main heatmap
            ),
            yaxis=dict(
                title='Fusion (Gene Pair & Breakpoint)',
                tickfont=dict(size=9),
                autorange='reversed'  # Top to bottom
            ),
            height=total_height,
            margin=dict(l=200, r=50, t=100, b=150),
            plot_bgcolor='white',
            font=dict(size=10)
        )
        
        # Save to HTML
        fig.write_html(
            output_html,
            include_plotlyjs='cdn',
            config={
                'displayModeBar': True,
                'displaylogo': False,
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': f'fusion_coverage_heatmap_{report_type.lower()}',
                    'height': total_height,
                    'width': 1400,
                    'scale': 2
                },
                'scrollZoom': True
            }
        )
        
        print(f"Heatmap visualization written to: {output_html}")
        
        # Print summary statistics
        print(f"\nSummary Statistics:")
        print(f"  Patients with 1 fusion: {np.sum(patient_fusion_counts == 1)}")
        print(f"  Patients with 2+ fusions: {np.sum(patient_fusion_counts >= 2)}")
        print(f"  Max fusions in one patient: {int(patient_fusion_counts.max())}")
        
        return True
        
    except Exception as e:
        print(f"Error creating heatmap: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description='Create fusion-patient coverage heatmap')
    parser.add_argument('--input', required=True,
                       help='Input fusion report TSV file')
    parser.add_argument('--output', required=True,
                       help='Output HTML file for heatmap visualization')
    parser.add_argument('--max-fusions', type=int, default=None,
                       help='Maximum number of fusions to display (default: all)')
    parser.add_argument('--report-type', choices=['Unthresholded', 'Thresholded'], 
                       default='Unthresholded',
                       help='Type of report being visualized (default: Unthresholded)')
    
    # New argument to optionally disable clustering
    parser.add_argument('--no-cluster', action='store_true',
                       help='Disable hierarchical clustering (will sort patients by fusion burden instead)')
    
    args = parser.parse_args()
    
    # Validate max_fusions
    if args.max_fusions is not None and args.max_fusions < 1:
        print("Error: --max-fusions must be at least 1")
        sys.exit(1)
    
    # Create heatmap
    success = create_fusion_heatmap(
        args.input,
        args.output,
        args.max_fusions,
        args.report_type,
        cluster=not args.no_cluster # Pass cluster boolean
    )
    
    if success:
        print("\nHeatmap visualization completed successfully!")
    else:
        print("\nHeatmap visualization skipped or failed.")
        sys.exit(0) # Don't fail the pipeline

if __name__ == "__main__":
    main()
