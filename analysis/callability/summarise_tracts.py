"""
create_snp_tract_bed.py - simple script to calculate inter-SNP tracts

intended to be run as part of callability.smk and not as a standalone script

from workflow:

input:
    expand('data/callability/tracts/{cross}.tracts.tsv', cross=CROSSES)
output:
    'data/callability/tract_distribution.tsv'

"""

import os
import csv
import itertools
from cyvcf2 import VCF
from tqdm import tqdm

def summarise_tracts(tract_files, tsv_output):
    with open(str(tsv_output), 'w') as f_out:
        fieldnames = ['cross', 'tract_length', 'count', 'sequence_length']
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for fname in tqdm(list(tract_files)):
            d = {}

            with open(fname, 'r') as f_in:
                reader = csv.DictReader(f_in, delimiter='\t')
                cross = os.path.basename(fname).rstrip('.tracts.tsv')

                for line in tqdm(reader, desc=cross):
                    length = int(line['length'])
                    if length not in d:
                        d[length] = 1
                    elif length in d:
                        d[length] += 1

                # fill in 'blanks'
                max_tract_length = max(d.keys())
                for tract_length in range(0, max_tract_length):
                    if tract_length not in d:
                        d[tract_length] = 0
                    else:
                        continue

                # write to file
                for tract_length in sorted(d.keys()):
                    writer.writerow({
                        'cross': cross,
                        'tract_length': tract_length,
                        'count': d[tract_length],
                        'sequence_length': int(tract_length * d[tract_length])
                        })
    
summarise_tracts(snakemake.input, snakemake.output)
