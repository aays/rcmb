"""
create_snp_tract_bed.py - simple script to calculate inter-SNP tracts

intended to be run as part of callability.smk and not as a standalone script

requires:
samtools

from workflow:

input:
    'data/alignments/parental_bam/{parent}.bam'
output:
    'data/callability/coverage/{parent}.cov.temp.tsv'

"""

import os
import csv
import subprocess
from tqdm import tqdm

lengths = {'chromosome_01': 8225636,
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
    'chromosome_17': 6954842}

def windowed_coverage(bam_fname, coverage_output, params, lengths):

    fieldnames = [
        '#rname', 'startpos', 'endpos', 'numreads', 'covbases', 
        'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
    cmd = 'samtools coverage -H -r {chrom}:{start}-{end} {fname}'
    bam_fname = str(bam_fname)

    out_fname = 'data/callability/coverage/' 
    out_fname += os.path.basename(bam_fname).rstrip('.bam')
    out_fname += '.cov.temp.tsv'
    with open(out_fname, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fieldnames)
        for chrom in tqdm(lengths):
            for window_start in tqdm(range(0, lengths[chrom], int(params.window_size / 2)), desc=chrom):
                window_end = window_start + params.window_size
                if window_end > lengths[chrom]:
                    window_end = lengths[chrom]
                cmd_split = cmd.format(
                        chrom=chrom,
                        start=window_start,
                        end=window_end,
                        fname=bam_fname).split(' ')
                proc = subprocess.run(
                    cmd_split,
                    check=True,
                    stdout=subprocess.PIPE)
                writer.writerow(
                    proc.stdout.decode('utf-8').rstrip('\n').split('\t')
                )

windowed_coverage(snakemake.input, snakemake.output, snakemake.params, lengths)
