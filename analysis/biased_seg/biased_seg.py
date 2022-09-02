"""
biased_seg.py - get counts of reads matching either ancestor across windows
"""

import time
import sys
import os
import re
import csv
import math
import subprocess
import argparse
import functools
import numpy as np
import pandas as pd
import pysam
import readcomb.classification as rc
from io import BytesIO
from multiprocessing import Pool
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
        description='get counts of reads matching either ancestor', 
        usage='python biased_seg.py [options]')

    parser.add_argument('-f', '--bam', required=True,
        type=str, help='Read name sorted sam/bam file containing reads of interest')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF with calls for cross of interest')
    parser.add_argument('-w', '--window_size', required=True, default=1000,
        type=int, help='Windowsize (default 1000)')
    parser.add_argument('-q', '--base_qual', required=False, default=0,
        type=int, help='Base qual filter (default 0)')
    parser.add_argument('-m', '--mapq', required=False, default=0,
        type=int, help='MAPQ filter (default 0)')
    parser.add_argument('-p', '--processes', required=False, default=1,
        type=int, help='Number of processes (default 1)')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    return parser

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

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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

def examine_pair(current_pair, lengths, args, header):
    """
    takes in a pair and returns which array to update and which window (e.g.
    parent1 array, chromosome_12, window idx 4 (e.g. window_size * 4)

    """
    record_1 = pysam.AlignedSegment.fromstring(current_pair[0], header)
    record_2 = pysam.AlignedSegment.fromstring(current_pair[1], header)
    pair = rc.Pair(record_1, record_2, args.vcf)
    pair.classify(masking=0, quality=args.base_qual)
    haps = list(set([sample for sample, _, _ in pair.condensed]))
    parents = re.search('[GB0-9]{4,5}x[0-9]{4}', args.bam).group(0).split('x')
    parents = ['CC' + parent for parent in parents if not parent.startswith('GB')]

    chrom = pair.rec_1.reference_name
    window = math.floor(int(pair.midpoint.split(':')[1]) / args.window_size)

    # ignore reads on contigs - no arrays to update
    if pair.rec_1.reference_name not in lengths.keys():
        return

    # mapq filter
    elif any([pair.rec_1.mapq < args.mapq, pair.rec_2.mapq < args.mapq]):
        return 'failed_mapq', chrom, window

    # remove anything with a phase change or no detectable variation
    elif (
        pair.call != 'no_phase_change' or
        not pair.condensed
    ):
        return 'removed', chrom, window

    elif pair.call == 'no_phase_change':
        assert len(haps) == 1
        hap = haps[0]
        if hap == parents[0]:
            return 'parent1', chrom, window
        elif hap == parents[1]:
            return 'parent2', chrom, window

    else:
        print(pair)
        raise Exception("this probably shouldn't happen")
        
# partial function with hardcoded args
# since args can't be passed to pool.imap otherwise
# ideally I'd use functools partial for this but apparently it's not pickleable???
# it also means I need to handle command line args globally, which I'm not happy about
parser = arg_parser()
args = parser.parse_args()
def examine_pair_wrapped(pair):
    header = pysam.AlignmentFile(args.bam, 'r').header
    return examine_pair(pair, lengths=lengths, args=args, header=header)

def biased_seg_parallel(args, lengths):
    """calculate biased seg stats using multiprocessing

    Parameters
    ----------
    args : argparse.ArgumentParser
        command line args
    lengths : dict
        dictionary containing chromosome lengths


    Returns
    -------
    None
    """

    print('[rcmb] creating processes')

    bam = pysam.AlignmentFile(args.bam, 'r')

    # get number of reads to track progress with
    stats = subprocess.check_output(['samtools', 'idxstats', args.bam])
    stats_table = pd.read_csv(
        BytesIO(stats), sep='\t',
        names=['chrom', 'length', 'map_reads', 'unmap_reads'])
    stats_table = stats_table[stats_table.chrom.isin(lengths.keys())]
    pair_count = int(stats_table['map_reads'].sum() / 2) # divide 2 to get pairs
    progress = tqdm(total=pair_count)

    # create arrays to store values in
    parent1_all = create_array_set(lengths, args.window_size)
    parent2_all = create_array_set(lengths, args.window_size)
    removed_all = create_array_set(lengths, args.window_size)
    failed_mapq_all = create_array_set(lengths, args.window_size)

    # https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
    # pull two reads at a time 
    def pairwise(iterable):
        reader = iter(iterable)
        for read1, read2 in zip(reader, reader):
            yield (read1.to_string(), read2.to_string())

    reader = pairwise(bam) # reads in two at once

    # increment arrays in parallel
    with Pool(processes=args.processes, maxtasksperchild=5000) as pool:
        for result in pool.imap(examine_pair_wrapped, reader):
            if not result:
                continue
            progress.update(n=1)
            call, chrom, window = result
            # f.write(str(result) + '\n')
            if call == 'parent1':
                parent1_all[chrom][window] += 1
            elif call == 'parent2':
                parent2_all[chrom][window] += 1
            elif call == 'removed':
                removed_all[chrom][window] += 1
            elif call == 'failed_mapq':
                failed_mapq_all[chrom][window] += 1
            else:
                raise ValueError(f'pair handled incorrectly - {call}, {chrom}, {window}')
                    
    # write to file
    with open(args.out, 'w') as f:
        fieldnames = [
            'chromosome', 'window_start', 'window_end', 'parent1_pairs', 'parent2_pairs', 
            'total_kept_pairs', 'parent1_perc', 'parent2_perc', 'pairs_skipped', 'failed_mapq']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        print('[rcmb] writing to file...')
        for chrom in tqdm(lengths):
            windows = range(0, lengths[chrom], args.window_size)
            for i, window_start in enumerate(windows):
                window_end = window_start + args.window_size \
                    if i < len(windows) - 1 else lengths[chrom]

                out_dict = {
                    'chromosome': chrom,
                    'window_start': window_start,
                    'window_end': window_end,
                    'parent1_pairs': int(parent1_all[chrom][i]),
                    'parent2_pairs': int(parent2_all[chrom][i]),
                    'pairs_skipped': int(removed_all[chrom][i]),
                    'failed_mapq': int(failed_mapq_all[chrom][i])
                    }

                out_dict['total_kept_pairs'] = int(parent1_all[chrom][i] + parent2_all[chrom][i])
                if not out_dict['total_kept_pairs'] == 0:
                    out_dict['parent1_perc'] = round(
                        out_dict['parent1_pairs'] / out_dict['total_kept_pairs'], 3)
                    out_dict['parent2_perc'] = round(
                        out_dict['parent2_pairs'] / out_dict['total_kept_pairs'], 3)
                else:

                    out_dict['parent1_perc'] = 'NA'
                    out_dict['parent2_perc'] = 'NA'

                writer.writerow(out_dict)

def main():
    biased_seg_parallel(args, lengths)

if __name__ == '__main__':
    main()

        

