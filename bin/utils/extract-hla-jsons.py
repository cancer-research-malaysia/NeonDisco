#!/usr/bin/env python3
import sys
import json

# grab the first command line argument
if len(sys.argv) != 3:
    print("Usage: python extract-hla-jsons.py <sample_name> <input_json>")
    sys.exit(1)
    
# assign the argument to a variable
sample_name = sys.argv[1]
input_json = sys.argv[2]

# read in json line from a file
with open(input_json, "r") as f:
    json_line = f.readline().strip()

# Parse the JSON
data = json.loads(json_line)

# Extract unique values for each key and format them
formatted_values = []

for key, values in data.items():
    # Get unique values by converting to a set
    unique_values = set(values)
    
    # Format each value with the "HLA-" prefix and process colons
    for value in unique_values:
        # Process the value to handle multiple colons
        parts = value.split(':')
        if len(parts) > 2:
            # Keep only the first two parts (remove second colon and trailing digits)
            processed_value = f"{parts[0]}:{parts[1]}"
        else:
            processed_value = value
            
        formatted_values.append(f"HLA-{processed_value}")

# Join the formatted values with commas
output_line = ",".join(formatted_values)

print(f"{sample_name}\t{output_line}")
