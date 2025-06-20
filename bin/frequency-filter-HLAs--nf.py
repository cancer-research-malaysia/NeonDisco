#!/usr/bin/env python3
"""
HLA Allotype Frequency Analysis
Extracts HLA allotypes with frequency ≥5% across the cohort for each HLA type (A, B, C)
"""

import sys
from collections import defaultdict, Counter
import argparse

def parse_hla_file(filename):
    """Parse HLA allotype file and return list of samples with their allotypes"""
    samples = []
    
    with open(filename, 'r') as f:
        # Skip header line
        next(f)
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            if len(parts) >= 2:
                sample_name = parts[0]
                allotypes = parts[1].split(',')
                # Clean up any whitespace
                allotypes = [allotype.strip() for allotype in allotypes if allotype.strip()]
                samples.append((sample_name, allotypes))
    
    return samples

def categorize_allotypes(samples):
    """Categorize allotypes by HLA type (A, B, C) and count frequencies"""
    hla_counts = {
        'HLA-A': Counter(),
        'HLA-B': Counter(), 
        'HLA-C': Counter()
    }
    
    total_samples = len(samples)
    
    for sample_name, allotypes in samples:
        for allotype in allotypes:
            # Determine HLA type
            if allotype.startswith('HLA-A*'):
                hla_counts['HLA-A'][allotype] += 1
            elif allotype.startswith('HLA-B*'):
                hla_counts['HLA-B'][allotype] += 1
            elif allotype.startswith('HLA-C*'):
                hla_counts['HLA-C'][allotype] += 1
    
    return hla_counts, total_samples

def filter_frequent_allotypes(hla_counts, total_samples, min_frequency=0.05):
    """Filter allotypes with frequency >= min_frequency (default 5%)"""
    frequent_allotypes = {}
    
    for hla_type, counts in hla_counts.items():
        frequent_allotypes[hla_type] = []
        
        for allotype, count in counts.items():
            frequency = count / total_samples
            if frequency >= min_frequency:
                frequent_allotypes[hla_type].append({
                    'allotype': allotype,
                    'count': count,
                    'frequency': frequency
                })
        
        # Sort by frequency (descending)
        frequent_allotypes[hla_type].sort(key=lambda x: x['frequency'], reverse=True)
    
    return frequent_allotypes

def print_results(frequent_allotypes, total_samples, min_frequency):
    """Print the results in a formatted way"""
    print(f"HLA Allotype Frequency Analysis")
    print(f"Total samples: {total_samples}")
    print(f"Minimum frequency threshold: {min_frequency*100:.1f}%")
    print("="*60)
    
    for hla_type in ['HLA-A', 'HLA-B', 'HLA-C']:
        print(f"\n{hla_type} allotypes with frequency ≥ {min_frequency*100:.1f}%:")
        print("-" * 50)
        
        if frequent_allotypes[hla_type]:
            for item in frequent_allotypes[hla_type]:
                print(f"{item['allotype']:<25} {item['count']:>3} samples ({item['frequency']*100:>5.1f}%)")
        else:
            print(f"No {hla_type} allotypes found with frequency ≥ {min_frequency*100:.1f}%")

def export_filtered_data(samples, frequent_allotypes, output_filename):
    """Export samples with only the frequent allotypes"""
    # Create set of frequent allotypes for quick lookup
    frequent_set = set()
    for hla_type, allotype_list in frequent_allotypes.items():
        for item in allotype_list:
            frequent_set.add(item['allotype'])
    
    with open(output_filename, 'w') as f:
        f.write("sampleName\thlaAllotypes\n")
        
        for sample_name, allotypes in samples:
            # Filter to only frequent allotypes
            filtered_allotypes = [allotype for allotype in allotypes if allotype in frequent_set]
            
            if filtered_allotypes:  # Only write if sample has frequent allotypes
                f.write(f"{sample_name}\t{','.join(filtered_allotypes)}\n")

def save_analysis_report(frequent_allotypes, total_samples, min_frequency, output_filename):
    """Save the analysis report to a file"""
    with open(output_filename, 'w') as f:
        f.write(f"HLA Allotype Frequency Analysis\n")
        f.write(f"Total samples: {total_samples}\n")
        f.write(f"Minimum frequency threshold: {min_frequency*100:.1f}%\n")
        f.write("="*60 + "\n")
        
        for hla_type in ['HLA-A', 'HLA-B', 'HLA-C']:
            f.write(f"\n{hla_type} allotypes with frequency ≥ {min_frequency*100:.1f}%:\n")
            f.write("-" * 50 + "\n")
            
            if frequent_allotypes[hla_type]:
                for item in frequent_allotypes[hla_type]:
                    f.write(f"{item['allotype']:<25} {item['count']:>3} samples ({item['frequency']*100:>5.1f}%)\n")
            else:
                f.write(f"No {hla_type} allotypes found with frequency ≥ {min_frequency*100:.1f}%\n")

def export_allotype_string(frequent_allotypes, output_filename):
    """Export frequent allotypes as a single comma-separated string"""
    # Collect all frequent allotypes in order: HLA-A, HLA-B, HLA-C
    all_frequent = []
    
    for hla_type in ['HLA-A', 'HLA-B', 'HLA-C']:
        for item in frequent_allotypes[hla_type]:
            all_frequent.append(item['allotype'])
    
    with open(output_filename, 'w') as f:
        f.write(','.join(all_frequent) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Filter HLA allotypes by frequency')
    parser.add_argument('input_file', help='Input HLA allotype file')
    parser.add_argument('--min-freq', type=float, default=0.05, 
                       help='Minimum frequency threshold (default: 0.05 = 5%%)')
    parser.add_argument('--output', '-o', help='Output file for analysis report')
    parser.add_argument('--string-output', help='Output file for comma-separated allotype string')
    
    args = parser.parse_args()
    
    # Generate default output filenames
    default_filtered_file = f"Frequency_filtered_HLA_allotypes.txt"
    default_report_file = f"Frequency_filtering_report.txt"
    default_string_file = f"HLA_allotypes_thresholded_cohortwide_string.txt"
    
    try:
        # Parse input file
        samples = parse_hla_file(args.input_file)
        
        if not samples:
            print("Error: No valid samples found in input file")
            sys.exit(1)
        
        # Analyze frequencies
        hla_counts, total_samples = categorize_allotypes(samples)
        frequent_allotypes = filter_frequent_allotypes(hla_counts, total_samples, args.min_freq)
        
        # Print results to console
        print_results(frequent_allotypes, total_samples, args.min_freq)
        
        # Always export filtered data (default behavior)
        export_filtered_data(samples, frequent_allotypes, default_filtered_file)
        print(f"\nFiltered sample data saved to: {default_filtered_file}")
        
        # Save analysis report (-o option or default)
        report_filename = args.output if args.output else default_report_file
        save_analysis_report(frequent_allotypes, total_samples, args.min_freq, report_filename)
        print(f"Analysis report saved to: {report_filename}")
        
        # Always export allotype string (default behavior, or custom filename)
        string_filename = args.string_output if args.string_output else default_string_file
        export_allotype_string(frequent_allotypes, string_filename)
        print(f"HLA allotype string saved to: {string_filename}")
        
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

