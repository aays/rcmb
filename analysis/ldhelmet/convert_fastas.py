"""
convert_fastas.py - simple script to convert per-sample fastas
into per-chr fastas with each containing a given chr across all samples

intended to be run as part of ldhelmet.smk, and not as a standalone script

from workflow: 

input:
    expand('data/ldhelmet/fasta/{sample}.fa', sample=SAMPLES)
output:
    expand('data/ldhelmet/fasta/chromosome_{n}.fa', n=list(range(1, 18)))

"""

import os
from Bio import SeqIO
from tqdm import tqdm

strains = {}
strains['NA1'] = ['CC3071', 'CC3086', 'GB119', 'CC1691', 'CC2935', 'CC3059', 'CC3062']
strains['NA2'] = ['CC2932', 'CC2344', 'CC2343', 'CC2931', 'CC2342', 'CC1952']
strains['ALL'] = strains['NA1'] + strains['NA2']

def convert_fastas(sample_fastas, out_fastas):
    for fname in tqdm(sample_fastas):
        sample_chrom = os.path.basename(fname).rstrip('.fa')
        sample, chrom = sample_chrom.split('.')
        population = [
            k for k in strains.keys() if k.startswith('NA')
            and sample in strains[k]][0]

        # should just be written to one of the two here
        with open(f'data/ldhelmet/fasta_chrom/{population}.{chrom}.fa', 'a') as f:
            rec = SeqIO.read(fname, 'fasta')
            f.write(f'>{sample}\n')
            f.write(f'{str(rec.seq)}\n')

        # always write to 'all' file
        with open(f'data/ldhelmet/fasta_chrom/ALL.{chrom}.fa', 'a') as f:
        
            reader = SeqIO.index(fname, 'fasta')
            f.write(f'>{sample}\n')
            f.write(f'{str(rec.seq)}\n')

convert_fastas(snakemake.input, snakemake.output)
    
    
