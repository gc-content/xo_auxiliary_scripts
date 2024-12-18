import argparse
import subprocess
import os
import tempfile
from collections import defaultdict
import sys


def run_bedtools_intersect(file1, file2, output_file, exclusive=False, wa=False, wb=False):
    """
    Use bedtools intersect to find overlaps (or exclusives) between two BED files.
    """
    bedtools_cmd = ["bedtools", "intersect", "-a", file1, "-b", file2]

    # Add -v for exclusive operation if requested
    if exclusive:
        bedtools_cmd.append("-v")

    # Add -wa or -wb flags if requested
    if wa:
        bedtools_cmd.append("-wa")
    if wb:
        bedtools_cmd.append("-wb")

    try:
        with open(output_file, "w") if output_file else sys.stdout as out:
            subprocess.run(bedtools_cmd, stdout=out, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools intersect: {e}", file=sys.stderr)


def process_files(file1, file2, output, exclusive, wa, wb):
    # Create temporary files grouped by the value in the fourth column
    temp_files1 = defaultdict(lambda: tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed'))
    temp_files2 = defaultdict(lambda: tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed'))

    # Read and write lines to the respective temporary files for file1
    with open(file1, 'r') as f1:
        for line in f1:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            value = parts[3]
            temp_files1[value].write(line)

    # Read and write lines to the respective temporary files for file2
    with open(file2, 'r') as f2:
        for line in f2:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            value = parts[3]
            temp_files2[value].write(line)

    # Close all the temp files to ensure data is written
    temp_files1_paths = {key: temp_file.name for key, temp_file in temp_files1.items()}
    temp_files2_paths = {key: temp_file.name for key, temp_file in temp_files2.items()}
    for temp_file in temp_files1.values():
        temp_file.close()
    for temp_file in temp_files2.values():
        temp_file.close()

    results = []

    # Perform intersection for each key in file1
    for value in temp_files1_paths:
        if value in temp_files2_paths:
            group_output = tempfile.NamedTemporaryFile(delete=False, suffix=".bed", mode="w").name
  
            # Run bedtools intersect or exclusive intersect
            run_bedtools_intersect(temp_files1_paths[value], temp_files2_paths[value], group_output, exclusive, wa, wb)

            # Read the result and append to final results
            with open(group_output, "r") as result_file:
                results.extend(result_file.readlines())

            # Clean up temporary output file
            os.remove(group_output)
        
        elif exclusive:
            # If exclusive and no corresponding group in file2, include all file1 records for this tag
            with open(temp_files1_paths[value], "r") as f1_group:
                results.extend(f1_group.readlines())

    # Clean up input temporary files
    for temp_file_path in temp_files1_paths.values():
        os.remove(temp_file_path)
    for temp_file_path in temp_files2_paths.values():
        os.remove(temp_file_path)

    # Write all results to the final output file, or print to stdout if output is None
    if output:
        with open(output, "w") as out_file:
            out_file.writelines(results)
    else:
        for line in results:
            print(line, end='')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Intersect genomic intervals between two files using bedtools, grouped by value.")
    parser.add_argument("-f1", "--file1", type=str, required=True, help="First input file")
    parser.add_argument("-f2", "--file2", type=str, required=True, help="Second input file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output file (default: stdout)")
    parser.add_argument("-v", "--exclusive", action="store_true", help="Run exclusive intersect (-v)")
    parser.add_argument("--wa", action="store_true", help="Report original A entries that have overlaps in B (-wa)")
    parser.add_argument("--wb", action="store_true", help="Report original B entries that overlap with A (-wb)")

    args = parser.parse_args()
    process_files(args.file1, args.file2, args.output, args.exclusive, args.wa, args.wb)
