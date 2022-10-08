"""
windowed_effective_sequence.py - get effective sequence for a given
cross but split into windows
"""

import os
import math
import csv
import argparse
import numpy as np
from tqdm import tqdm
import readcomb.classification as rc

def args():
    parser = argparse.ArgumentParser(
        description='get effective sequence in windows', 
        usage='python windowed_effective_sequence.py [options]')

    parser.add_argument('-b', '--bam', required=True,
        type=str, help='BAM file containing all reads considered for phase changes')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='Corresponding VCF file')
    parser.add_argument('-w', '--window_size', required=True,
        type=int, help='Window size (default 1000)')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.bam, args.vcf, args.window_size, args.out

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
    'chromosome_17': 6954842,
    'cpDNA': 205535,
    'mtDNA': 15789 # including organelles for mt leakage
    }

def create_array_set(lengths, window_size):
    """generate generic set of arrays for provided chromosomes/contigs
    given a window size, will create length-L numpy arrays for each provided chromosome
    where L is chrom_length / window_size. the individual elements of the array
    represent windows, and windows can be recapitulated via element_index * window_size
    Parameters
    ----------
    lengths : dict
        dict of chromosome/contig/scaffold lengths
    window_size : int
        desired window size to split chrs/contigs/scaffolds into
    Returns
    -------
    arrays : dict
        dict of numpy arrays, with keys corresponding to chrs/contigs/scaffolds
    """
    arrays = {}
    for chrom in lengths.keys():
        array_len = len(range(0, lengths[chrom], window_size)) + 1 # add 1 for window end
        arrays[chrom] = np.zeros(array_len)
    return arrays

def parse_pair(pair):
    """get effective sequence + midpoint of given pair

    Parameters
    ----------
    pair : classification.Pair
        read pair of interest

    Returns
    -------
    chrom : str
        pair chromosome
    midpoint : int
        pair midpoint (halfway point)
    effective_length : int
        total effective sequence of pair
    """
    chrom = pair.rec_1.reference_name
    start = pair.rec_1.reference_start
    end = pair.rec_2.reference_start + len(pair.segment_2)
    effective_length = end - start
    midpoint = start + int(effective_length / 2) # will round down if needed
    return chrom, midpoint, effective_length

def denominator_counts(bam, vcf, window_size, out, lengths):
    """get effective sequence per window

    Parameters
    ----------
    bam : str
        path to readcomb prepped BAM file for cross
    vcf : str
        path to readcomb prepped VCF file for cross
    window_size : int
        window size to use
    out : str
        file to write to

    Returns
    -------
    None
    """
    reader = rc.pairs_creation(bam, vcf)
    effective_seq_arrays = create_array_set(lengths, window_size)
    read_count_arrays = create_array_set(lengths, window_size)
    cross = os.path.basename(bam).rstrip('.sorted.bam')

    for pair in tqdm(reader, desc=cross):
        if pair.rec_1.reference_name not in lengths: # ignore non chroms
            continue
        chrom, midpoint, effective_seq = parse_pair(pair)
        window = math.floor(midpoint / window_size)
        effective_seq_arrays[chrom][window] += effective_seq
        read_count_arrays[chrom][window] += 1

    with open(out, 'w') as f:
        fieldnames = [
            'cross', 'chromosome', 'window_start', 'window_end',
            'read_count', 'effective_sequence']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for chrom in tqdm(lengths):
            windows = range(0, lengths[chrom], window_size)
            for i, window_start in enumerate(windows):
                window_end = window_start + window_size \
                    if i < len(windows) - 1 else lengths[chrom]

                out_dict = {
                    'cross': cross,
                    'chromosome': chrom,
                    'window_start': window_start,
                    'window_end': window_end,
                    'read_count': int(read_count_arrays[chrom][i]),
                    'effective_sequence': int(effective_seq_arrays[chrom][i])
                    }
                writer.writerow(out_dict)

def main():
    bam, vcf, window_size, out = args()
    denominator_counts(bam, vcf, window_size, out, lengths)

if __name__ == '__main__':
    main()

        

