#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd

def parse_fusion_transcripts(fusion_dict):
    """
    Parse fusion transcript IDs from a dictionary of {sampleID: [transcripts]}
    
    Args:
    fusion_dict (dict): Dictionary with sampleID as keys and lists of fusion transcript IDs as values
    
    Returns:
    dict: Nested dictionary with sampleID as keys and lists of command dictionaries as values
    """
    # Initialize results dictionary
    parsed_commands = {}
    
    # Process each sample
    for sample_id, transcripts in fusion_dict.items():
        # Create nested dictionary for this sample
        parsed_commands[sample_id] = []
        
        # Process each transcript for this sample
        for transcript in transcripts:
            try:
                # Split the transcript into gene pair and genomic locations
                parts = transcript.split('__')
                if len(parts) != 2:
                    print(f"Warning: Malformed transcript ID: {transcript}. Skipping.")
                    continue
                    
                gene_pair, genomic_locations = parts
                
                # Split gene pair
                gene_parts = gene_pair.split('::')
                if len(gene_parts) != 2:
                    print(f"Warning: Malformed gene pair: {gene_pair}. Skipping.")
                    continue
                    
                gene1, gene2 = gene_parts
                
                # Split genomic locations
                location_parts = genomic_locations.split('-')
                if len(location_parts) != 2:
                    print(f"Warning: Malformed genomic locations: {genomic_locations}. Skipping.")
                    continue
                    
                chrom_pos1, chrom_pos2 = location_parts
                
                # Extract chromosomes and positions
                cp1_parts = chrom_pos1.split(':')
                cp2_parts = chrom_pos2.split(':')
                
                if len(cp1_parts) != 2 or len(cp2_parts) != 2:
                    print(f"Warning: Malformed chromosome:position format. Skipping.")
                    continue
                
                chrom1, pos1 = cp1_parts
                chrom2, pos2 = cp2_parts
                
                # Construct command dictionary with sampleID included
                command_dict = {
                    'sampleID': sample_id,
                    'gene5': gene1,
                    'gene3': gene2,
                    'junction5': pos1,
                    'junction3': pos2,
                    'output': f"{sample_id}_{gene1}-{gene2}"  # Include sampleID in output name
                }
                
                parsed_commands[sample_id].append(command_dict)
            except Exception as e:
                print(f"Error parsing transcript {transcript}: {e}")
    
    return parsed_commands


def generate_agfusion_commands_for_sample(parsed_commands, sample_id, db_path='/tmp/agfusion.homo_sapiens.111.db', output_dir="./agfusion-dirs"):
    """
    Generate all AGFusion bash commands for a specific sample ID
    
    Args:
    parsed_commands (dict): Nested dictionary from parse_fusion_transcripts
    sample_id: ID of the sample to generate commands for
    db_path (str): Path to the AGFusion database
    output_dir (str, optional): Base directory for output files
    
    Returns:
    list: List of bash commands for the specified sample
    """
    if sample_id not in parsed_commands:
        return []
        
    commands = []
    for cmd_dict in parsed_commands[sample_id]:
        output_path = f"{output_dir}/{cmd_dict['output']}"
            
        command = f"agfusion annotate \\\n" \
                 f"  -g5 {cmd_dict['gene5']} \\\n" \
                 f"  -g3 {cmd_dict['gene3']} \\\n" \
                 f"  -j5 {cmd_dict['junction5']} \\\n" \
                 f"  -j3 {cmd_dict['junction3']} \\\n" \
                 f"  -db {db_path} \\\n" \
                 f"  -o {output_path} \\\n" \
                 f" --middlestar \\\n" \
                 f"  --noncanonical"
        commands.append(command)
    
    return commands


def generate_all_commands(parsed_results, db_path='/tmp/agfusion.homo_sapiens.111.db', output_dir="./agfusion-dirs"):
    """Generate commands for all samples"""
    all_commands = {}
    for sample_id in parsed_results:
        commands = generate_agfusion_commands_for_sample(
            parsed_results, sample_id, db_path, output_dir)
        all_commands[sample_id] = commands
    
    return all_commands


def write_bash_script(output_file, all_commands, log_dir, output_prefix):
    """
    Write commands to a bash script with error handling
    
    Args:
    output_file (str): Path to output bash script
    all_commands (dict): Dictionary of {sample_id: [commands]}
    log_dir (str, optional): Directory to store logs
    output_prefix (str): Base directory to store agfusion output folders
    """
    with open(output_file, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated AGFusion commands with error handling\n\n")
        
        # Create log directory
        if log_dir:
            f.write(f"# Create log directory\n")
            f.write(f"mkdir -p {log_dir}\n\n")
            f.write(f"log_dir=\"{log_dir}\"\n\n")
        else:
            log_dir = "agfusion-logs"
            f.write(f"# Create log directory\n")
            f.write(f"mkdir -p {log_dir}\n\n")
            f.write(f"log_dir=\"{log_dir}\"\n\n")

        # Create output directory
        f.write(f"# Create output base directory\n")
        f.write(f"mkdir -p {output_prefix}\n\n")
        
        # Initialize counters
        f.write("# Initialize counters\n")
        f.write("TOTAL=0\n")
        f.write("SUCCESS=0\n")
        f.write("FAILED=0\n\n")
        
        # Helper function to run commands with error handling
        f.write("# Function to run command and handle errors\n")
        f.write("run_command() {\n")
        f.write("  local cmd=\"$1\"\n")
        f.write("  local desc=\"$2\"\n")
        f.write("  local log_file=\"$3\"\n")
        f.write("  \n")
        f.write("  echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Running: $desc\"\n")
        f.write("  echo \"$ $cmd\" >> \"$log_file\"\n")
        f.write("  \n")
        f.write("  TOTAL=$((TOTAL+1))\n")
        f.write("  if eval \"$cmd\" >> \"$log_file\" 2>&1; then\n")
        f.write("    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $desc\" | tee -a \"$log_file\"\n")
        f.write("    SUCCESS=$((SUCCESS+1))\n")
        f.write("    return 0\n")
        f.write("  else\n")
        f.write("    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] FAILED: $desc\" | tee -a \"$log_file\"\n")
        f.write("    echo \"Check $log_file for details\"\n")
        f.write("    FAILED=$((FAILED+1))\n")
        f.write("    return 1\n")
        f.write("  fi\n")
        f.write("}\n\n")
        
        # Write commands for each sample
        for sample_id, commands in all_commands.items():
            f.write(f"# Commands for Sample ID: {sample_id}\n")
            f.write(f"echo \"Processing sample {sample_id} with {len(commands)} fusions\"\n")
            
            for i, cmd in enumerate(commands, 1):
                # Extract gene pair for description
                gene5 = cmd.split("-g5 ")[1].split(" ")[0]
                gene3 = cmd.split("-g3 ")[1].split(" ")[0]
                
                # Create a unique log file for each command
                log_file = f"{log_dir}/{sample_id}_{gene5}-{gene3}.log"
                
                # Write a descriptive label for this command
                desc = f"sample {sample_id}, fusion {i}/{len(commands)}: {gene5}-{gene3}"
                
                # Use the run_command function to execute with error handling
                f.write(f"run_command \"{cmd}\" \"{desc}\" \"{log_file}\"\n\n")
        
        # Summary at the end
        f.write("# Print summary\n")
        f.write("echo \"\"\n")
        f.write("echo \"=== EXECUTION SUMMARY ===\"\n")
        f.write("echo \"Total commands: $TOTAL\"\n")
        f.write("echo \"Successful: $SUCCESS\"\n")
        f.write("echo \"Failed: $FAILED\"\n")
        f.write("echo \"Logs directory: $(realpath $log_dir)\"\n")
        f.write("echo \"=========================\"\n")
        f.write("\n")
        f.write("# Return overall success/failure\n")
        f.write("if [ \"$FAILED\" -eq 0 ]; then\n")
        f.write("  echo \"All commands completed successfully\"\n")
        f.write("  exit 0\n")
        f.write("else\n")
        f.write("  echo \"Some commands failed, check logs for details\"\n")
        f.write("  # Return 0 to avoid stopping the pipeline, but warn about failures\n")
        f.write("  # Change to 'exit 1' if you want the pipeline to stop on any failures\n")
        f.write("  exit 0\n")
        f.write("fi\n")
    
    # Make the script executable
    os.chmod(output_file, 0o755)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process fusion data and generate AGFusion commands')
    parser.add_argument('--input', '-i', required=True, help='Input TSV file with fusion data')
    parser.add_argument('--sample', '-s', help='Specific sample ID to process (optional)')
    parser.add_argument('--output-dir', '-o', default='./agfusion-dirs', help='Output base directory for AGFusion results')
    parser.add_argument('--db-path', '-d', default='/tmp/agfusion.homo_sapiens.111.db', 
                        help='Path to AGFusion database')
    parser.add_argument('--output-script', '-c', action='store_true',
                        help='Write commands to bash script file')
    parser.add_argument('--log-dir', '-l', default='agfusion-logs',
                        help='Directory to store command logs (default: agfusion-logs)')
    
    args = parser.parse_args()
    
    try:
        # Read input file
        print(f"Reading fusion data from {args.input}")
        fusion_data = pd.read_csv(args.input, sep='\t')
        
        # Group by sampleID and get unique fusionTranscriptIDs
        fusion_dict = fusion_data.groupby('sampleID')['fusionTranscriptID'].unique().to_dict()
        
        # Convert numpy arrays to lists
        for key in fusion_dict:
            fusion_dict[key] = fusion_dict[key].tolist()
        
        # Parse transcripts
        parsed_results = parse_fusion_transcripts(fusion_dict)
        
        # Generate commands
        if args.sample:
            # For a specific sample
            if args.sample not in parsed_results:
                print(f"Error: Sample ID {args.sample} not found in the data")
                sys.exit(1)
                
            commands = generate_agfusion_commands_for_sample(
                parsed_results, args.sample, args.db_path, args.output_dir)
            
            sample_commands = {args.sample: commands}
            
            if args.output_script:
                # Write to bash script with error handling
                write_bash_script(args.output_script, sample_commands, args.log_dir, args.output_dir)
                print(f"Bash script written to {args.output_script}")
            else:
                # Print to stdout
                print(f"\nCommands for Sample ID {args.sample}:")
                for cmd in commands:
                    print(cmd)
                    print()  # Add blank line between commands
        else:
            # For all samples
            all_commands = generate_all_commands(
                parsed_results, args.db_path, args.output_dir)
            
            # Output commands
            if args.output_script:
                # Write to bash script with error handling
                bash_script_name ='agfusion-cmd.sh'
                write_bash_script(bash_script_name, all_commands, args.log_dir, args.output_dir)
                print(f"Bash script written to {bash_script_name}")
            else:
                # Print to stdout
                for sample_id, cmds in all_commands.items():
                    print(f"\nCommands for Sample ID {sample_id}:")
                    for cmd in cmds:
                        print(cmd)
                        print()  # Add blank line between commands
                    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
