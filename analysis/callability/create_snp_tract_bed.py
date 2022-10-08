"""
create_snp_tract_bed.py - simple script to calculate inter-SNP tracts

intended to be run as part of callability.smk and not as a standalone script

from workflow:

input:
    'data/genotyping/vcf_filtered/{cross}.vcf.gz'
output:
    'data/callability/tracts/{cross}.tracts.tsv'
"""

import os
import csv
import itertools
from cyvcf2 import VCF
from tqdm import tqdm

# from itertools documentation
# since I don't have py3.10
def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def create_snp_tract_bed(vcf_input, tsv_output, wildcards):
    reader = VCF(str(vcf_input))
    pairwise_reader = pairwise(reader)
    with open(str(tsv_output), 'w') as f:
        fieldnames = [
            'cross', 'chrom', 'snp_1', 'snp_2', 'length']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        print(vcf_input)
        for snp_1, snp_2 in tqdm(pairwise_reader):
            if snp_1.CHROM != snp_2.CHROM: # edges of chrom
                continue
            writer.writerow({
                'cross': wildcards.cross,
                'chrom': snp_1.CHROM,
                'snp_1': snp_1.POS,
                'snp_2': snp_2.POS,
                'length': snp_2.POS - snp_1.POS
                })
    
create_snp_tract_bed(snakemake.input, snakemake.output, snakemake.wildcards)
