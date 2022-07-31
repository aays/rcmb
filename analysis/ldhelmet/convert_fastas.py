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

def convert_fastas(sample_fastas, out_fastas):
    for i in tqdm(range(1, 18), desc='chroms'):
        n = f'0{i}' if i < 10 else f'{i}'
        chrom = f'chromosome_{n}'
        
        for fname in tqdm(sample_fastas, desc='samples'):
            sample = os.path.basename(fname).rstrip('.fa')
            
            with open(f'data/ldhelmet/fasta_chrom/chromosome_{i}.fa', 'a') as f:
                
                reader = SeqIO.index(fname, 'fasta')
                f.write(f'>{sample}\n')
                f.write(f'{str(reader[chrom].seq)}\n')

convert_fastas(snakemake.input, snakemake.output)
    
    
