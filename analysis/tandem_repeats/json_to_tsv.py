"""
json_to_tsv.py - simple script to convert json output from ULTRA into
per-sample-chr tsv files

intended to be run as part of tandem_repeats_calc.smk and not as a standalone
script

from workflow:

input:
    'data/tandem_repeats/json/{sample}.{chr}.json'
output:
    'data/tandem_repeats/tsv/{sample}.{chr}.tsv'
"""

import os
import csv
import json
from tqdm import tqdm

def json_to_tsv(json_input, csv_output):
    for json_file in json_input:
        outfile_name = os.path.basename(json_file).replace('.json', '.tsv')
        with open(json_file, 'r') as f_in:
            with open('data/tandem_repeats/tsv/' + outfile_name, 'w') as f_out:
                reader = json.load(f_in)
                fieldnames = list(reader['Repeats'][0].keys())
                writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
                writer.writeheader()

                for record in tqdm(reader['Repeats']):
                    writer.writerow(record)

json_to_tsv(snakemake.input, snakemake.output)

