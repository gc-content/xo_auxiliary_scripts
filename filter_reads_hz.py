import pandas as pd
import pysam
import argparse
import sys
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# Helper function to handle defaultdict of lists
def defaultdict_of_list():
    return defaultdict(list)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter reads based on overlapping heterozygous variants.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file with variants")
    parser.add_argument("-o", "--output", help="Output file for filtered BAM")
    parser.add_argument("-t", "--threshold", type=float, default=0.5, help="Proportion threshold for heterozygous variants (default: 0.5)")
    parser.add_argument("-p", "--procs", type=int, default=cpu_count(), help="Number of processes (default: all CPUs)")
    parser.add_argument("--stats_output", required=True, help="Output file for the table of read proportions")
    parser.add_argument("--min_variants", type=int, default=3, help="Minimum number of variants to consider a read (default: 3)")
    return parser.parse_args()

# Custom VCF parser function that extracts sample names and organizes variants
def parse_vcf(vcf_file):
    variants = defaultdict(defaultdict_of_list)
    sample_names = []

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                # Extract sample names from the header line
                sample_names = line.strip().split('\t')[9:]
            elif line.startswith('#'):
                continue
            else:
                # Parsing each line of the VCF
                chrom, pos, _, ref, alt, _, _, _, fmt, *samples = line.strip().split('\t')
                pos = int(pos)

                # Process genotype (GT) for each sample
                for idx, sample in enumerate(samples):
                    gt = sample.split(':')[0]  # Get the genotype field
                    if gt == './.':  # Skip missing data
                        continue

                    hp1, hp2 = gt.split('/')

                    var = {
                        'is_het': (hp1 == '0' and hp2 == '1') or (hp1 == '1' and hp2 == '0')  # Check if heterozygous
                    }

                    # Store the variant in the structure for the correct chromosome and sample
                    sample_name = sample_names[idx]
                    variants[chrom][sample_name].append((pos, var))

    return variants

# Function to process each chromosome in the BAM file
def process_chromosome(chrom, variants, bam_file, threshold):
    results = []
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Dictionary to track the last processed variant index for each sample (LN tag)
    last_processed_index = defaultdict(int)

    # Process each read for the current chromosome
    for read in bam.fetch(chrom):
        ln_tag = dict(read.tags).get('LN')
        if ln_tag is None or ln_tag not in variants[chrom]:
            continue  # Skip reads without an LN tag or no matching variants

        read_start = read.reference_start
        read_end = read.reference_end

        # Get the variants for the specific line (LN tag) in the current chromosome
        sample_variants = variants[chrom][ln_tag]

        total_variants = 0
        het_variants = 0

        # Get the last processed index for the current sample
        idx = last_processed_index[ln_tag]

        # Skip variants that are before the read start
        while idx < len(sample_variants) and sample_variants[idx][0] < read_start:
            idx += 1

        # Continue processing variants that overlap the read
        for variant_idx in range(idx, len(sample_variants)):
            var_pos, var_info = sample_variants[variant_idx]
            if var_pos > read_end:
                break  # Variants beyond the read end, stop processing

            total_variants += 1
            if var_info['is_het']:
                het_variants += 1

        # Calculate the proportion of heterozygous variants
        if total_variants > 0:
            proportion_het = het_variants / total_variants
        else:
            proportion_het = 0

        # Update the last processed index for the current sample
        last_processed_index[ln_tag] = idx

        results.append({
            "read_name": read.query_name,
            "chromosome": chrom,
            "read_start": read_start,
            "read_end": read_end,
            "ln_tag": ln_tag,
            "total_variants": total_variants,
            "het_variants": het_variants,
            "proportion_het": proportion_het
        })

    bam.close()
    return results

# Function for parallel processing across chromosomes
def parallel_process(chromosomes, variants, bam_file, threshold, num_procs):
    with Pool(processes=num_procs) as pool:
        results = pool.starmap(process_chromosome, [(chrom, variants, bam_file, threshold) for chrom in chromosomes])
    return [item for sublist in results for item in sublist]

# Main function
def main(bam_file, vcf_file, output_file, threshold, stats_output, num_procs, min_variants=3):
    # Step 1: Parse the VCF file and organize variants by chromosome and sample
    variants_dict = parse_vcf(vcf_file)

    # Step 2: Get all chromosomes from the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    chromosomes = bam.references
    bam.close()

    # Step 3: Process BAM file with parallelization
    results = parallel_process(chromosomes, variants_dict, bam_file, threshold, num_procs)

    # Step 4: Output results (proportion table for plotting)
    results_df = pd.DataFrame(results)
    results_df.to_csv(stats_output, index=False, sep="\t")

    # Step 5: Filter BAM based on threshold
    if output_file:
        bam_in = pysam.AlignmentFile(bam_file, "rb")
        bam_out = pysam.AlignmentFile(output_file, "wb", header=bam_in.header)

        # Collect reads that exceed the threshold or don't meet minimum variant requirements (reads to remove)
        reads_to_remove = {row['read_name'] for _, row in results_df.iterrows()
                           if row['proportion_het'] > threshold and row['total_variants'] >= min_variants}

        for read in bam_in.fetch():
            if read.query_name not in reads_to_remove:
                bam_out.write(read)

        bam_in.close()
        bam_out.close()

if __name__ == "__main__":
    args = parse_arguments()
    main(args.bam, args.vcf, args.output, args.threshold, args.stats_output, args.procs, args.min_variants)
