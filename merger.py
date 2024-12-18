import argparse
from collections import Counter, defaultdict
import sys


def parse_data(data):
    """
    Parse the input data into a list of tuples.
    Each tuple is of the form (chr, start, end, value).
    """
    parsed = []
    for line in data.strip().split("\n"):
        parts = line.split()
        chr = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        value = int(parts[3])
        parsed.append((chr, start, end, value))
    return parsed


def merge_intervals(intervals, use_broader):
    """
    Merge overlapping intervals based on the most common start/end,
    or the broader/narrower interval in case of a tie.
    """
    if not intervals:
        return []

    # Sort intervals by start position
    intervals.sort(key=lambda x: x[1])

    merged = []
    current_start, current_end, current_chr, current_value = intervals[0][1], intervals[0][2], intervals[0][0], intervals[0][3]
    start_list = [current_start]
    end_list = [current_end]
    count = 1

    for i in range(1, len(intervals)):
        chr, start, end, value = intervals[i]

        if chr == current_chr and start <= current_end:  # Overlapping intervals
            start_list.append(start)
            end_list.append(end)
            current_end = max(current_end, end)
            count += 1
        else:
            # Resolve start and end based on majority vote or chosen strategy (broader/narrower)
            resolved_start = resolve_coordinate(start_list, conservative=not use_broader)
            resolved_end = resolve_coordinate(end_list, conservative=use_broader)
            merged.append((current_chr, resolved_start, resolved_end, current_value, count))

            # Reset for the next set of intervals
            current_chr, current_value = chr, value
            current_start, current_end = start, end
            start_list = [start]
            end_list = [end]
            count = 1

    # Handle the last group
    resolved_start = resolve_coordinate(start_list, conservative=not use_broader)
    resolved_end = resolve_coordinate(end_list, conservative=use_broader)
    merged.append((current_chr, resolved_start, resolved_end, current_value, count))

    return merged


def resolve_coordinate(coord_list, conservative):
    """
    Resolve the majority vote for start or end coordinates; if there's a tie, use the conservative approach.
    If conservative is True, return the lowest start; if False, return the highest end.
    """
    counter = Counter(coord_list)
    most_common = counter.most_common()

    # Check if there's a tie
    if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
        return max(coord_list) if conservative else min(coord_list)

    # Return the most common value
    return most_common[0][0]


def main(input_stream, output_stream, use_broader):
    # Read data from input
    data = input_stream.read()

    # Parse the data
    parsed_data = parse_data(data)

    # Group data by chromosome and value (column 4)
    grouped_data = defaultdict(list)
    for entry in parsed_data:
        chr, start, end, value = entry
        grouped_data[(chr, value)].append(entry)

    # Process each group to merge intervals
    results = []
    for (chr, value), intervals in grouped_data.items():
        merged_intervals = merge_intervals(intervals, use_broader)
        results.extend(merged_intervals)

    # Write results to output
    for result in results:
        output_stream.write("\t".join(map(str, result)) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge overlapping genomic intervals by chromosome and value.")
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), default=sys.stdin,
                        help="Input file (default: stdin)")
    parser.add_argument("-o", "--output", type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file (default: stdout)")
    parser.add_argument("-b", "--broader", action="store_true",
                        help="Use the broader interval in case of tie (default: narrower)")

    args = parser.parse_args()
    main(args.input, args.output, args.broader)
