"""
ldhelmet workflow - generates per-sample fastas from previously
generated GVCFs, converts them to per-chr fastas, and then runs
LDhelmet on them

requires:
bcftools, gatk4, ldhelmet

usage:
snakemake -pr -s analysis/ldhelmet/ldhelmet.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

SAMPLES = list(set([
    os.path.basename(fname.rstrip('.bam'))
    for fname in glob('data/alignments/parental_bam/*')
    if fname.endswith('.bam')]
))

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/ldhelmet/fasta')
mkdir('data/ldhelmet/fasta_chrom')

# --- rules

rule all:
    input:
        expand('data/ldhelmet/fasta/{sample}.fa', sample=SAMPLES),
        expand('data/ldhelmet/fasta_chrom/chromosome_{i}.fa', i=list(range(1, 18)))

rule create_sample_fasta:
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa',
        gvcf = 'data/genotyping/gvcf_sample/{sample}.g.vcf.gz'
    output:
        'data/ldhelmet/fasta/{sample}.fa'
    log:
        'data/ldhelmet/fasta/log.txt'
    threads:
        4
    shell:
        'time bcftools consensus --fasta-ref {input.ref} --output {output} '
        '--sample {wildcards.sample} {input.gvcf} 2>> {log}'

rule create_per_chr_fastas:
    input:
        expand('data/ldhelmet/fasta/{sample}.fa', sample=SAMPLES)
    output:
        expand('data/ldhelmet/fasta_chrom/chromosome_{n}.fa', n=list(range(1, 18)))
    script:
        'convert_fastas.py'


