"""
tsv_to_bed.py - simple script to convert per-chr tsvs into per sample files

intended to be run as part of tandem_repeats_calc.smk and not as a standalone
script

from workflow:

input:
    'data/tandem_repeats/tsv/{sample}.{chr}.json'
output:
    'data/tandem_repeats/samples/{sample}.tsv'
"""

import os
import re
import csv
import json
from tqdm import tqdm

SAMPLES = [
    'CC1691', 'CC1952', 'CC2342', 'CC2343', 'CC2344', 'CC2931',
    'CC2932', 'CC2935', 'CC3059', 'CC3062', 'CC3071', 'CC3086', 'GB119']

def tsv_to_bed(tsv_input, sample_output):

    fieldnames = [
        '#chrom', 'start', 'end', 'length', 'period', 'score', 'log2_pval', 
        'substitutions', 'insertions', 'deletions', 'consensus',
        'consensus_length', 'sequence']

    for sample in SAMPLES:

        with open(f'data/tandem_repeats/samples/{sample}.tsv', 'w') as f_out:
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()

            sample_files = sorted([fname for fname in tsv_input if sample in fname])
            for sample_file in tqdm(sample_files):
                with open(sample_file, 'r') as f_in:
                    reader = csv.DictReader(f_in, delimiter='\t')
                    for line in reader:
                        out_dict = {
                        '#chrom': re.search('chromosome_[0-9]{2}', sample_file).group(0),
                        'start': line['Start'],
                        'end': int(line['Start']) + int(line['Length']),
                        'length': int(line['Length']),
                        'period': line['Period'],
                        'score': line['Score'],
                        'log2_pval': line['Log2 Pval'],
                        'substitutions': line['Substitutions'],
                        'insertions': line['Insertions'],
                        'deletions': line['Deletions'],
                        'consensus': line['Consensus'],
                        'consensus_length': len(line['Consensus']),
                        'sequence': line['Sequence']
                        }
                        writer.writerow(out_dict)


tsv_to_bed(snakemake.input, snakemake.output)

