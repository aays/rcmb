"""
snp_density.py - simple script to get windowed SNP counts

intended to be run as part of rate.smk and not as a standalone script

from workflow:

input:
    'data/genotyping/vcf_filtered/{cross}.vcf.gz'
output:
    'data/rate/snp_density/{cross}.snps.2kb.tsv'
"""

import os
import csv
from cyvcf2 import VCF
from tqdm import tqdm

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

def windowed_snp_count(vcf_input, tsv_output, wildcards, params, lengths):
    reader = VCF(str(vcf_input))
    window_size = int(params.window_size)
    with open(str(tsv_output), 'w') as f:
        fieldnames = [
            'cross', 'chromosome', 'window_start', 'window_end', 'snp_count', 'snp_density']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for chrom in tqdm(lengths.keys(), desc=wildcards.cross):
            windows = range(0, lengths[chrom], window_size)
            for i, window_start in tqdm(enumerate(windows), desc=chrom):
                window_end = window_start + window_size \
                    if i < len(windows) - 1 else lengths[chrom]

                snp_count = len([snp for snp in reader(f'{chrom}:{window_start}-{window_end}')])

                out_dict = {
                    'cross': wildcards.cross,
                    'chromosome': chrom,
                    'window_start': window_start,
                    'window_end': window_end,
                    'snp_count': snp_count,
                    'snp_density': round(snp_count / (window_end - window_start), 3)
                }
                writer.writerow(out_dict)

windowed_snp_count(snakemake.input, snakemake.output, snakemake.wildcards, snakemake.params, lengths)
