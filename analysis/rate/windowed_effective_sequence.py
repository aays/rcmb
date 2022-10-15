"""
windowed_effective_sequence.py - get effective sequence for a given
cross but split into windows
"""

import os
import math
import csv
import argparse
import pysam
from tqdm import tqdm
import readcomb.classification as rc

def args():
    parser = argparse.ArgumentParser(
        description='get effective sequence in windows', 
        usage='python windowed_effective_sequence.py [options]')

    parser.add_argument('-b', '--bam', required=True,
        type=str, help='BAM file containing all reads considered for phase changes')
    parser.add_argument('-s', '--snp_read_counts', required=True,
        type=str, help='SNP tract read count file')
    parser.add_argument('-w', '--window_size', required=True,
        type=int, help='Window size (default 1000)')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.bam, args.snp_read_counts, args.window_size, args.out

lengths = {
    'chromosome_01': 8225636,
    'chromosome_02': 8655884,
    'chromosome_03': 9286894,
    'chromosome_04': 4130073,
    'chromosome_05': 3682160,
    'chromosome_06': 8913359,
    'chromosome_07': 6492107,
    'chromosome_08': 4526983,
    'chromosome_09': 6807148,
    'chromosome_10': 6800247,
    'chromosome_11': 4479522,
    'chromosome_12': 9952739,
    'chromosome_13': 5281438,
    'chromosome_14': 4217303,
    'chromosome_15': 5870643,
    'chromosome_16': 8042475,
    'chromosome_17': 6954842
    }

def get_eff_bp(bam_reader, chrom, start, end):
    eff_bp = 0
    for p_column in bam_reader.pileup(chrom, start, end, truncate=True):
        for p_read in p_column.pileups:
            eff_bp += 1
    return eff_bp

def denominator_counts(bam, snp_read_counts, window_size, out, lengths):

    chrom_windows = {}
    for chrom in lengths:
        window_range = 0, lengths[chrom]
        chrom_windows[chrom] = list(range(0, lengths[chrom], window_size))
        chrom_windows[chrom].append(lengths[chrom])

    with open(snp_read_counts, 'r') as f_in:
        with open(out, 'w') as f_out:
            reader = csv.DictReader(f_in, delimiter='\t')
            fieldnames = [
                'cross', 'chrom', 'window_start', 'window_end', 'eff_bp']
            cross = os.path.basename(bam.rstrip('sorted.bam'))
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()

            reader = pysam.TabixFile(snp_read_counts)
            bam_reader = pysam.AlignmentFile(bam)

            for chrom in tqdm(chrom_windows):
                windows = chrom_windows[chrom]
                for i, window_start in tqdm(enumerate(windows), total=len(windows)):
                    window_end = window_start + window_size \
                        if i < len(windows) - 1 else lengths[chrom]

                    tracts = [
                        r for r in reader.fetch(chrom, window_start, window_end)]
                    coords = []

                    # get tracts
                    for tract in tracts:
                        _, start, end, _, pairs_span, _ = tract.split('\t')
                        start, end, pairs_span = int(start), int(end), int(pairs_span)
                        if pairs_span > 0:
                            coords.append([start, end])
                    
                    # handle
                    if not coords:
                        out_dict = {
                            'cross': cross,
                            'chrom': chrom,
                            'window_start': window_start,
                            'window_end': window_end,
                            'eff_bp': 0
                            }
                        writer.writerow(out_dict)
                    
                    window_eff_bp = 0
                    for snp_1, snp_2 in coords:
                        if snp_1 < window_start:
                            window_eff_bp += get_eff_bp(
                                bam_reader, chrom, window_start, snp_2)
                        elif snp_2 > window_end:
                            window_eff_bp += get_eff_bp(
                                bam_reader, chrom, snp_1, window_end)
                        else:
                            window_eff_bp += get_eff_bp(
                                bam_reader, chrom, snp_1, snp_2)
                    out_dict = {
                        'cross': cross,
                        'chrom': chrom,
                        'window_start': window_start,
                        'window_end': window_end,
                        'eff_bp': window_eff_bp
                        }
                    writer.writerow(out_dict)


def main():
    bam, snp_read_counts, window_size, out = args()
    denominator_counts(bam, snp_read_counts, window_size, out, lengths)

if __name__ == '__main__':
    main()

        

