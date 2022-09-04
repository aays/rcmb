"""
tandem repeats calc workflow - will generate
tabixed bed files of repeats

requires:
singularity, ULTRA

usage:
snakemake -pr -s analysis/biased_seg/tandem_repeats_calc.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

SAMPLES = list(set([
    os.path.basename(fname.rstrip('.bam'))
    for fname in glob('data/alignments/parental_bam/*')
    if fname.endswith('.bam')]
))

CHROMS = [f'chromosome_0{i}' for i in range(1, 10)]
CHROMS.extend([f'chromosome_{i}' for i in range(10, 18)])

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/tandem_repeats/json')
mkdir('data/tandem_repeats/tsv')
mkdir('data/tandem_repeats/samples')

# --- rules

rule all:
    input:
        expand('data/tandem_repeats/json/{sample}.{chrom}.json', 
            sample=SAMPLES, chrom=CHROMS),
        expand('data/tandem_repeats/tsv/{sample}.{chrom}.tsv',
            sample=SAMPLES, chrom=CHROMS),
        expand('data/tandem_repeats/samples/{sample}.sorted.tsv.gz', 
            sample=SAMPLES),
        expand('data/tandem_repeats/samples/{sample}.sorted.tsv.gz.tbi',
            sample=SAMPLES)

rule ultra_run:
    input:
        'data/ldhelmet/fasta/{sample}.{chrom}.fa'
    output:
        'data/tandem_repeats/json/{sample}.{chrom}.json'
    threads:
        8
    params:
        at_content = '0.36',
        score_threshold = '9',
        min_mapq = '1',
        abs_path = '/research/projects/chlamydomonas/genomewide_recombination/rcmb'
    shell:
        'time singularity exec --bind '
        '{params.abs_path}/data/ldhelmet/fasta_chrom/:'
        '{params.abs_path}/analysis/tandem_repeats/mnt '
        'analysis/tandem_repeats/ultra.sif '
        'ultra -f {output} -at {params.at_content} '
        '-n {threads} -s {params.score_threshold} {input}'

rule json_to_tsv:
    input:
        expand('data/tandem_repeats/json/{sample}.{chrom}.json', sample=SAMPLES, chrom=CHROMS)
    output:
        expand('data/tandem_repeats/tsv/{sample}.{chrom}.tsv', sample=SAMPLES, chrom=CHROMS)
    script:
        'json_to_tsv.py'

rule tsv_to_bed:
    input:
        expand('data/tandem_repeats/tsv/{sample}.{chrom}.tsv', sample=SAMPLES, chrom=CHROMS)
    output:
        expand('data/tandem_repeats/samples/{sample}.tsv', sample=SAMPLES)
    script:
        'tsv_to_bed.py'

rule sort_bed:
    input:
        temp('data/tandem_repeats/samples/{sample}.tsv')
    output:
        temp('data/tandem_repeats/samples/{sample}.sorted.tsv')
    shell:
        'sort -k1,1 -k2n,3n {input} > {output}'

rule bgzip_bed:
    input:
        'data/tandem_repeats/samples/{sample}.sorted.tsv'
    output:
        'data/tandem_repeats/samples/{sample}.sorted.tsv.gz'
    shell:
        'bgzip {input}'

rule tabix_bed:
    input:
        'data/tandem_repeats/samples/{sample}.sorted.tsv.gz'
    output:
        'data/tandem_repeats/samples/{sample}.sorted.tsv.gz.tbi'
    shell:
        'tabix -p bed -b2 -e3 {input}'
    
