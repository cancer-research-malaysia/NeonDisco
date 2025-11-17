#!/usr/bin/env python3
"""
Script to create interactive heatmap of patient-fusion coverage
Generates both interactive HTML (Plotly) and static PDF (PyComplexHeatmap)
"""
import sys
import argparse
import polars as pl
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

def create_plotly_heatmap(df, matrix_final, fusion_labels_final, patients_final, 
                          patient_fusion_counts_final, output_html, report_type, 
                          total_height, all_patients):
    """Create interactive Plotly heatmap"""
    # Create hover text
    hover_text = []
    for row in df.iter_rows(named=True):
        hover_row = []
        patients_with_fusion = set(row['SamplesWithFusion'].split(','))
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
        z=matrix_final,
        x=patients_final,
        y=fusion_labels_final,
        hovertext=hover_text,
        hoverinfo='text',
        colorscale=[[0, '#f0f0f0'], [1, '#d62728']],
        showscale=False,
        xgap=1,
        ygap=1
    ))
    
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
        ),
        yaxis=dict(
            title='Fusion (Gene Pair & Breakpoint)',
            tickfont=dict(size=9),
            autorange='reversed'
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
    
    print(f"Interactive HTML heatmap written to: {output_html}")


def create_pycomplexheatmap_pdf(df_final, matrix_final, fusion_labels_final, patients_final,
                                patient_fusion_counts_final, output_pdf, report_type,
                                all_patients):
    """Create publication-quality PDF heatmap using PyComplexHeatmap"""
    try:
        from PyComplexHeatmap import ClusterMapPlotter, HeatmapAnnotation, anno_simple
        
        # Try to use fastcluster for better performance
        try:
            import fastcluster
            print("Using fastcluster for improved clustering performance")
        except ImportError:
            print("Note: fastcluster not installed, using scipy (slower). Install with: pip install fastcluster")
        
        # Create simplified fusion labels - use full breakpoint for uniqueness
        fusion_labels_simple = []
        for row in df_final.iter_rows(named=True):
            # Use full breakpoint ID to ensure uniqueness
            label = f"{row['fusionGenePair']}|{row['breakpointID']}"
            fusion_labels_simple.append(label)
        
        # Create pandas DataFrame (required by PyComplexHeatmap)
        df_heatmap = pd.DataFrame(
            matrix_final,
            index=fusion_labels_simple,
            columns=patients_final
        )
        
        # Create column annotation for patient fusion burden
        df_col_anno = pd.DataFrame({
            'Fusion Count': patient_fusion_counts_final
        }, index=patients_final)
        
        # Create top annotation showing fusion burden
        col_ha = HeatmapAnnotation(
            **{'Fusion\nCount': anno_simple(
                df_col_anno['Fusion Count'],
                cmap='Blues',
                legend=True,
                height=5
            )},
            axis=1,
            verbose=0,
            label_side='left',
            label_kws={'horizontalalignment': 'right'}
        )
        
        # Calculate figure size dynamically
        n_fusions = len(df_final)
        n_patients = len(patients_final)
        fig_width = max(10, min(24, n_patients * 0.35))
        fig_height = max(8, min(32, n_fusions * 0.28))
        
        # Create the clustermap
        plt.figure(figsize=(fig_width, fig_height))
        cm = ClusterMapPlotter(
            data=df_heatmap,
            top_annotation=col_ha,
            row_cluster=True,
            col_cluster=True,
            row_cluster_method='average',
            row_cluster_metric='jaccard',
            col_cluster_method='average',
            col_cluster_metric='jaccard',
            row_dendrogram=True,
            col_dendrogram=True,
            show_rownames=True,
            show_colnames=True,
            row_names_side='right',
            col_names_side='bottom',
            cmap='Reds',  # White for absence (0), red for presence (1)
            label='Presence',
            legend=True,  # Enable heatmap legend
            legend_side='right',  # Position legends on the right
            legend_hpad=10,  # Horizontal padding from heatmap (mm)
            legend_vpad=5,  # Vertical padding from top (mm)
            row_split_gap=0.5,
            col_split_gap=0.5,
            xticklabels_kws={'labelrotation': -90, 'labelsize': 7},
            yticklabels_kws={'labelsize': 6},
            rasterized=True,
            vmin=0,
            vmax=1,
            verbose=0
        )
        
        # Add title
        plt.suptitle(
            f'{report_type} Recurrent Fusion Coverage Matrix\n'
            f'{len(all_patients)} patients × {len(df_final)} fusions',
            fontsize=14,
            y=0.98
        )
        
        # Save to PDF
        plt.savefig(output_pdf, bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"Static PDF heatmap written to: {output_pdf}")
        return True
        
    except ImportError:
        print("\nWarning: PyComplexHeatmap not installed. Skipping PDF generation.")
        print("To install: pip install PyComplexHeatmap")
        return False
    except Exception as e:
        print(f"\nWarning: Could not create PDF heatmap: {e}")
        import traceback
        traceback.print_exc()
        return False


def create_fusion_heatmap(report_tsv, output_html, output_pdf=None, max_fusions=None, 
                         report_type="Unthresholded", cluster=True):
    """
    Create interactive heatmap showing which patients have which fusions
    
    Args:
        report_tsv: Path to fusion report TSV file
        output_html: Output HTML file path
        output_pdf: Output PDF file path (optional)
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
        
        # Create fusion labels (gene pair + breakpoint for HTML display)
        fusion_labels = []
        for row in df.iter_rows(named=True):
            if row:
                # Use full breakpoint ID for uniqueness
                label = f"{row['fusionGenePair']}<br><sub>{row['breakpointID']}</sub>"
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
                if patient in patient_to_idx:
                    patient_idx = patient_to_idx[patient]
                    matrix[fusion_idx, patient_idx] = 1
        
        # Calculate patient fusion burden
        patient_fusion_counts = matrix.sum(axis=0)

        # Initialize default order variables
        fusion_order = list(range(len(df)))
        patient_order = list(range(len(all_patients)))
        
        # Perform hierarchical clustering if requested
        if cluster and len(df) >= 2 and len(all_patients) >= 2:
            print("Performing hierarchical clustering...")
            
            # Cluster rows (fusions)
            fusion_distances = pdist(matrix, metric='jaccard')
            fusion_linkage = linkage(fusion_distances, method='average')
            from scipy.cluster.hierarchy import dendrogram
            fusion_dendrogram = dendrogram(fusion_linkage, no_plot=True)
            fusion_order = fusion_dendrogram['leaves']
            
            # Cluster columns (patients)
            patient_distances = pdist(matrix.T, metric='jaccard')
            patient_linkage = linkage(patient_distances, method='average')
            patient_dendrogram = dendrogram(patient_linkage, no_plot=True)
            patient_order = patient_dendrogram['leaves']

            print(f"Clustering complete: {len(df)} fusions × {len(all_patients)} patients")
        else:
            print("Sorting by fusion burden (not clustering)...")
            patient_order = np.argsort(-patient_fusion_counts)
        
        # Apply final order
        matrix_final = matrix[fusion_order, :][:, patient_order]
        fusion_labels_final = [fusion_labels[i] for i in fusion_order]
        patients_final = [all_patients[i] for i in patient_order]
        patient_fusion_counts_final = patient_fusion_counts[patient_order]
        df_final = df[fusion_order]
        
        # Determine height
        base_height = 400
        height_per_fusion = 25
        total_height = max(base_height, base_height + len(df) * height_per_fusion)
        
        # Create Plotly interactive heatmap
        create_plotly_heatmap(
            df_final, matrix_final, fusion_labels_final, patients_final,
            patient_fusion_counts_final, output_html, report_type,
            total_height, all_patients
        )
        
        # Create PyComplexHeatmap PDF if requested
        if output_pdf:
            create_pycomplexheatmap_pdf(
                df_final, matrix_final, fusion_labels_final, patients_final,
                patient_fusion_counts_final, output_pdf, report_type, all_patients
            )
        
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
    parser = argparse.ArgumentParser(
        description='Create fusion-patient coverage heatmap (HTML + PDF)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate both HTML and PDF
  %(prog)s --input fusions.tsv --output heatmap.html --pdf heatmap.pdf
  
  # HTML only
  %(prog)s --input fusions.tsv --output heatmap.html
  
  # Limit to top 50 fusions
  %(prog)s --input fusions.tsv --output heatmap.html --max-fusions 50
        """
    )
    parser.add_argument('--input', required=True,
                       help='Input fusion report TSV file')
    parser.add_argument('--output', required=True,
                       help='Output HTML file for interactive heatmap')
    parser.add_argument('--pdf', default=None,
                       help='Output PDF file for static heatmap (requires PyComplexHeatmap)')
    parser.add_argument('--max-fusions', type=int, default=None,
                       help='Maximum number of fusions to display (default: all)')
    parser.add_argument('--report-type', choices=['Unthresholded', 'Thresholded'], 
                       default='Unthresholded',
                       help='Type of report being visualized (default: Unthresholded)')
    parser.add_argument('--no-cluster', action='store_true',
                       help='Disable hierarchical clustering')
    
    args = parser.parse_args()
    
    # Validate max_fusions
    if args.max_fusions is not None and args.max_fusions < 1:
        print("Error: --max-fusions must be at least 1")
        sys.exit(1)
    
    # Create heatmap
    success = create_fusion_heatmap(
        args.input,
        args.output,
        args.pdf,
        args.max_fusions,
        args.report_type,
        cluster=not args.no_cluster
    )
    
    if success:
        print("\nHeatmap visualization completed successfully!")
    else:
        print("\nHeatmap visualization skipped or failed.")
        sys.exit(0)


if __name__ == "__main__":
    main()

    