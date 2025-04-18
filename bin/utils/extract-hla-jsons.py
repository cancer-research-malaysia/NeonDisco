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
    
    # Format each value with the "HLA-" prefix
    for value in unique_values:
        formatted_values.append(f"HLA-{value}")

# Join the formatted values with commas
output_line = ",".join(formatted_values)

# Write to a TSV file
with open(f"{sample_name}-HLA-types-reformatted.tsv", "w") as f:
    f.write(output_line)

print(f"Output: {output_line}")
