import vcf
import numpy as np
import pandas as pd
import argparse
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter variants based on allele frequency and sample count criteria.")
    parser.add_argument("vcf_file", type=str, help="Input VCF file path.")
    parser.add_argument("-o", "--output", type=str,
                        help="Output file to save the filtered BED file. Default is stdout.")
    parser.add_argument("-s", "--min_samples", type=int, default=10,
                        help="Minimum number of samples in which a variant must be present to be included. Default is 10.")
    parser.add_argument("-m", "--metric", type=str, choices=['mean', 'median'], default='mean',
                        help="Statistic to calculate for allele frequency: 'mean' or 'median'. Default is 'mean'.")
    return parser.parse_args()


def calculate_mean_allele_frequency(freqs):
    valid_freqs = [f for f in freqs if f is not None]
    return np.mean(valid_freqs) if valid_freqs else 0


def calculate_median_allele_frequency(freqs):
    valid_freqs = [f for f in freqs if f is not None]
    return np.median(valid_freqs) if valid_freqs else 0


def filter_variants(vcf_file, output, min_samples, metric):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    bed_records = []

    # Choose the appropriate function based on the metric
    if metric == 'mean':
        calculate_af = calculate_mean_allele_frequency
        metric_label = "Mean_AF"
    else:
        calculate_af = calculate_median_allele_frequency
        metric_label = "Median_AF"

    for record in vcf_reader:
        # Count samples with a non-None AF value
        sample_count = sum(1 for call in record.samples if call.data.AF is not None)

        # Check if variant meets the minimum samples criterion
        if sample_count >= min_samples:
            freqs = [call.data.AF for call in record.samples if call.data.AF is not None]
            af_value = calculate_af(freqs)
            bed_records.append((record.CHROM, record.POS - 1, record.POS, af_value))

    df_bed = pd.DataFrame(bed_records, columns=["Chromosome", "Start", "End", metric_label])

    if output:
        df_bed.to_csv(output, sep="\t", header=False, index=False)
    else:
        df_bed.to_csv(sys.stdout, sep="\t", header=False, index=False)


if __name__ == "__main__":
    args = parse_arguments()
    filter_variants(args.vcf_file, args.output, args.min_samples, args.metric)