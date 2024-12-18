import pandas as pd
import pysam
import argparse
import sys
from multiprocessing import Pool, cpu_count

# Function to process each chromosome
def process_chromosome(chrom, bed_df, bam_file, s=0.9, t=0.1):
    results = []
    # Filter BED file for current chromosome
    chrom_variants = bed_df[bed_df['chrom'] == chrom].reset_index(drop=True)
    variant_idx = 0  # Start index for variants in the chromosome

    # Open BAM file for reading
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Iterate through each read in the BAM file for the given chromosome
    for read in bam.fetch(chrom):
        read_start = read.reference_start
        read_end = read.reference_end

        # Find overlapping variants using the locking approach
        overlapping_variants = []
        
        # Move the index up to the first variant that overlaps the read or beyond
        while variant_idx < len(chrom_variants) and chrom_variants.at[variant_idx, 'end'] < read_start:
            variant_idx += 1
        
        # Collect all variants that overlap the read
        temp_idx = variant_idx
        while temp_idx < len(chrom_variants) and chrom_variants.at[temp_idx, 'start'] <= read_end:
            overlapping_variants.append(chrom_variants.iloc[temp_idx])
            temp_idx += 1

        # If there are overlapping variants, calculate metrics
        if overlapping_variants:
            overlapping_variants = pd.DataFrame(overlapping_variants)
            below_threshold_count = (overlapping_variants['score'] < s).sum()
            proportion_below_threshold = below_threshold_count / len(overlapping_variants)
            average_score = overlapping_variants['score'].mean()

            # Append results
            results.append({
                "read_name": read.query_name,
                "chromosome": chrom,
                "read_start": read_start,
                "read_end": read_end,
                "proportion_below_threshold": proportion_below_threshold,
                "average_score": average_score
            })
        else:
            results.append({
                "read_name": read.query_name,
                "chromosome": chrom,
                "read_start": read_start,
                "read_end": read_end,
                "proportion_below_threshold": "",
                "average_score": None
            })

    bam.close()
    return results

# Multiprocessing wrapper
def parallel_process(chromosomes, bed_df, bam_file, s, t, num_procs):
    with Pool(processes=num_procs) as pool:
        results = pool.starmap(process_chromosome, [(chrom, bed_df, bam_file, s, t) for chrom in chromosomes])
    return [item for sublist in results for item in sublist]  # Flatten the list of results

# Main function
def main(bam_file, bed_file, output_file, bam_out, s, t, num_procs):
    # Load the BED file (assuming columns: chrom, start, end, score)
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "score"])

    # Get all chromosomes from BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    chromosomes = bam.references
    bam.close()

    # Perform parallel processing to get results
    results = parallel_process(chromosomes, bed_df, bam_file, s, t, num_procs)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Output results either to a file or to stdout
    if output_file:
        results_df.to_csv(output_file, index=False, sep="\t")
    else:
        results_df.to_csv(sys.stdout, index=False, sep="\t")

    # Open BAM file for reading and filtered BAM file for writing
    if bam_out:
        # Prepare for writing the filtered BAM
        bam_in = pysam.AlignmentFile(bam_file, "rb")
        filtered_bam = pysam.AlignmentFile(bam_out, "wb", header=bam_in.header)
        
        # Create a lookup dictionary for the reads to be written
        reads_to_write = {result['read_name']: result for result in results if result['proportion_below_threshold'] != "" and result['proportion_below_threshold'] <= t}

        # Iterate through reads and write if they are in the filtered set
        for read in bam_in.fetch():
            if read.query_name in reads_to_write:
                filtered_bam.write(read)
        
        bam_in.close()
        filtered_bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate variant score statistics for reads in BAM file.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-v", "--bed", required=True, help="Input BED file with variants and scores")
    parser.add_argument("-o", "--output", help="Output file for stats (default: stdout)")
    parser.add_argument("--bam_out", help="Output BAM file for filtered reads")
    parser.add_argument("-s", "--score_threshold", type=float, default=0.9, help="Score threshold for filtering (default: 0.9)")
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="Proportion threshold for reads to keep (default: 0.1)")
    parser.add_argument("-p", "--procs", type=int, default=cpu_count(), help="Number of processes (default: all CPUs)")

    args = parser.parse_args()

    # Run the main function
    main(args.bam, args.bed, args.output, args.bam_out, args.score_threshold, args.threshold, args.procs)

