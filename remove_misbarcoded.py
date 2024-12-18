import pysam
import argparse
from collections import defaultdict, deque
from multiprocessing import Pool

def main():
    parser = argparse.ArgumentParser(description="Filter misbarcoded reads from BAM file based on haplotype tagging.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("--window", type=int, default=20000, help="Window size (in bp) around the read to check for haplotyped reads (default: 20,000 bp)")
    parser.add_argument("--n_reads", type=int, default=5, help="Number of haplotyped reads required in the window to apply filtering (default: 5)")
    parser.add_argument("--num_processes", type=int, default=1, help="Number of processes to use (default: 1)")
    
    args = parser.parse_args()

    # Process BAM file by chromosomes with multiprocessing
    process_bam_by_chromosomes(args.input_bam, args.window, args.n_reads, args.num_processes)

def process_bam_by_chromosomes(input_bam, window_size, n_reads, num_processes):

    # Open BAM file and get list of chromosomes
    bam = pysam.AlignmentFile(input_bam, "rb")
    chromosomes = bam.references
    bam.close()

    # Use multiprocessing to process each chromosome separately
    with Pool(processes=num_processes) as pool:
        results = []
        for chrom in chromosomes:
            results.append(pool.apply_async(process_bam_chromosome, (input_bam, chrom, window_size, n_reads)))

        # Wait for all processes to finish
        for result in results:
            result.wait()


def process_bam_chromosome(input_bam, chrom, window_size, n_reads):
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Store reads by individual and haplotype in a deque, sorted by end position
    window_reads = defaultdict(lambda: deque())

    # Set to track processed reads
    processed_reads = defaultdict(lambda: deque())

    # Iterate over all reads in the current chromosome
    for read in bam.fetch(chrom):
        # Skip secondary and supplementary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        # Skip reads that don't have the required tags (LN and HP)
        if not read.has_tag("LN") or not read.has_tag("HP"):
            continue

        individual = read.get_tag("LN")
        haplotype = read.get_tag("HP")
        read_start = read.reference_start
        read_end = read.reference_end  # End of the current read

        # Check the p structure for upstream reads
        process_previous_reads(window_reads[individual], read_start, window_size, n_reads, read, processed_reads)

        # Add current read to p (sorted by end)
        add_read_to_window(window_reads[individual], (read_end, read))

    # Process remaining reads after loop
    flush_remaining_reads(window_reads, n_reads, chrom, processed_reads)

    bam.close()

def process_previous_reads(window_reads, current_start, window_size, n_reads, current_read, processed_reads):
    """Process previous reads whose window ends before the current read's start."""
    # Use a temporary list to hold reads that are in the window
    in_window_reads = []

    individual = current_read.get_tag("LN")  # Get individual from the current read

    # Iterate through the existing window reads without removing them
    for end_pos, read_to_process in window_reads:
        if end_pos + window_size < current_start:
            in_window_reads.append(read_to_process)

    # Process the reads that are in the window of the current read
    for read_to_process in in_window_reads:
        if read_to_process not in processed_reads[individual]:  # Check against the individual-specific deque
            process_window(read_to_process, window_reads, window_size, n_reads, processed_reads)


def process_window(read, window_reads, window_size, n_reads, processed_reads):
    """Process reads within the window of a given read."""
    read_start = read.reference_start
    haplotype = read.get_tag("HP")
    individual = read.get_tag("LN")

    # Count haplotyped reads in the window
    haplotype_counts = defaultdict(int)

    # Pop reads that are outside the window from the left
    while window_reads and window_reads[0][1].reference_end < read_start - window_size:
        window_reads.popleft()  # Maintain the window correctly

    # Count the remaining reads only if they haven't been processed yet
    for end_pos, stored_read in window_reads:
        haplotype_counts[stored_read.get_tag("HP")] += 1

    # If haplotype reads are sufficient, check for misbarcoding
    total_haplotype_reads = sum(haplotype_counts.values())
    if total_haplotype_reads >= n_reads and haplotype_counts[haplotype] == 1:
        # Print the misbarcoded read information in TSV format
        print(f"{read.query_name}\t{read.reference_name}\t{read.reference_start}\t{read.reference_end}\t{individual}")

    # Mark the current read as processed for the correct individual
    processed_reads[individual].append(read)  # Use the individual-specific deque

    # Pop outdated reads from processed_reads for the current individual
    while processed_reads[individual] and processed_reads[individual][0].reference_end < read_start - window_size:
        processed_reads[individual].popleft()  # Maintain the deque correctly for this individual

def add_read_to_window(window_reads, new_read):
    """Insert the new read into the sorted window by end coordinate."""
    read_end, read = new_read

    # Find the position to insert the new read to keep the list sorted by end position
    if not window_reads or window_reads[-1][0] <= read_end:
        window_reads.append(new_read)
    else:
        for i in range(len(window_reads)-1, -1, -1):
            if window_reads[i][0] <= read_end:
                window_reads.insert(i+1, new_read)
                break
        else:
            window_reads.appendleft(new_read)

def flush_remaining_reads(window_reads, n_reads, chrom, processed_reads):
    """Flush any remaining reads after processing."""
    for individual, reads in window_reads.items():
        for _, read in reads:
            # Process the remaining reads as needed
            process_window(read, reads, 20000, n_reads, processed_reads)

if __name__ == "__main__":
    main()

