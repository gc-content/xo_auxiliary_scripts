import pandas as pd
from intervaltree import Interval, IntervalTree
import argparse
import sys

# Argument parser setup
parser = argparse.ArgumentParser(description="Intersect or find non-overlapping intervals between two files.")
parser.add_argument("file1", help="Path to the first file (e.g., file1.txt) with intervals.")
parser.add_argument("file2", help="Path to the second file (e.g., file2.txt) with intervals and detailed format.")
parser.add_argument("--output", help="Path for the output file. Defaults to stdout.", default=None)
parser.add_argument("-v", "--invert", action="store_true", help="Output non-intersecting rows instead of intersecting ones.")
args = parser.parse_args()

# Load data from file1 and file2 with all columns for file2
file1 = pd.read_csv(args.file1, sep="\t", header=None, names=["chromosome", "start", "end", "individual_id", "score"])
file2 = pd.read_csv(args.file2, sep="\t", header=None, 
                    names=["ReadID", "Line", "Chromosome", "StartPos", "IsRecombinant", "PreviousReadID",
                           "HaplotypeBlocks", "VariantString", "HaplotypeString", "Breakpoints", "SwitchScore", "LocalSwitchScore"])

# Convert necessary columns to integer type
file1["start"] = file1["start"].astype(int)
file1["end"] = file1["end"].astype(int)
file1["individual_id"] = file1["individual_id"].astype(int)
file2["Line"] = file2["Line"].astype(int)

# Parsing intervals in the 10th column (interval_data) in file2
def parse_intervals(interval_str):
    intervals = []
    for part in interval_str.split(';'):
        if part:
            range_part = part.split(":")[0]  # Extract `start-end` part
            start, end = map(int, range_part.split("-"))
            intervals.append((start, end))
    return intervals

file2["intervals"] = file2["Breakpoints"].apply(parse_intervals)

# Create a dictionary of IntervalTrees for each chromosome and individual in File 1
trees = {}
for _, row in file1.iterrows():
    chrom = row["chromosome"]
    ind_id = row["individual_id"]
    key = (chrom, ind_id)
    
    if key not in trees:
        trees[key] = IntervalTree()
    # Storing the entire row directly in the interval tree for easier access
    trees[key][row["start"]:row["end"]] = row
    
# Replace missing or NA fields with 'NA' to ensure they are printed properly
file2["PreviousReadID"] = file2["PreviousReadID"].fillna("NA")  # Replace NaN with 'NA'

# Find intersections or non-intersections based on `-v` flag
results = []
for _, row in file2.iterrows():
    chrom = row["Chromosome"]
    ind_id = row["Line"]
    key = (chrom, ind_id)

    has_overlap = False
    if key in trees:
        for start, end in row["intervals"]:
            intervals = trees[key].overlap(start, end)
            if intervals:
                has_overlap = True
                break

    # Add row to results based on `-v` flag logic
    if (has_overlap and not args.invert) or (not has_overlap and args.invert):
        # Only add columns present in the original `file2`
        row_data = row.drop("intervals")  # Drop the extra intervals column for output
        results.append(row_data)

# Convert results to DataFrame in the exact format of file2
intersection_df = pd.DataFrame(results, columns=file2.columns)

# Output to stdout or file based on args.output
if args.output:
    intersection_df.to_csv(args.output, sep="\t", index=False)
else:
    intersection_df.to_csv(sys.stdout, sep="\t", index=False)

