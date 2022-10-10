"""
summarise_chroms.py - summarise chrom callables

intended to be run as part of callability.smk and not as a standalone script

from workflow:

input:
    expand('data/callability/tracts/{cross}.read_counts.tsv', cross=CROSSES)
output:
    'data/callability/callables/chrom_callables.tsv'

"""

import os
import csv
import itertools
from cyvcf2 import VCF
from tqdm import tqdm

def summarise_tracts(tract_files, tsv_output):
    with open(str(tsv_output), 'w') as f_out:
        fieldnames = [
            'cross', 'chrom', 'total_sequence', 'callable_sequence', 'snp_count']
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for fname in tqdm(list(tract_files)):
            d = {}

            with open(fname, 'r') as f_in:
                reader = csv.DictReader(f_in, delimiter='\t')
                cross = os.path.basename(fname).rstrip('.read_counts.tsv')

                for line in tqdm(reader, desc=cross):
                    chrom = line['chrom']
                    if chrom not in d:
                        d[chrom] = {}
                        d[chrom]['total_sequence'] = 0
                        d[chrom]['snp_count'] = 0
                        d[chrom]['callable_sequence'] = 0
                    d[chrom]['total_sequence'] += int(line['length'])
                    d[chrom]['snp_count'] += 1
                    if int(line['pairs_span']) > 0:
                        d[chrom]['callable_sequence'] += int(line['length'])

                # write to file
                for chrom in sorted(d.keys()):
                    writer.writerow({
                        'cross': cross,
                        'chrom': chrom,
                        'total_sequence': d[chrom]['total_sequence'],
                        'callable_sequence': d[chrom]['callable_sequence'],
                        'snp_count': d[chrom]['snp_count']
                        })
    
summarise_tracts(snakemake.input, snakemake.output)
