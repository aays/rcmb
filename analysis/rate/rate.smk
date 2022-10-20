"""
rate calc workflow

requires:
readcomb

usage:
snakemake -pr -s analysis/rate/rate.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

with open('data/genotyping/samples.txt', 'r') as f:
    CROSSES = [line.split(' ')[0] for line in f]
    CROSSES.remove('3071x2931')
    CROSSES.remove('3071x3062')
    CROSS_DICT = {}
    for line in CROSSES:
        mt_plus, mt_minus = line.split('x')
        if mt_plus not in CROSS_DICT:
            CROSS_DICT[mt_plus] = [mt_minus]
        else:
            CROSS_DICT[mt_plus].append(mt_minus)

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/rate/denominators/')
mkdir('data/rate/denominators/2kb')
mkdir('data/rate/snp_density/')
mkdir('data/callability/callables/annotations')

# --- rules

rule all:
    input:
        expand('data/rate/denominators/{cross}.tsv', cross=CROSSES),
        expand('data/rate/snp_density/{cross}.snps.2kb.tsv', cross=CROSSES),
        expand('data/callability/callables/annotations/{cross}.tsv', cross=CROSSES)

rule windowed_denominator_calc:
    input:
        bam = 'data/alignments/bam_prepped/sorted/{cross}.sorted.bam',
        tracts = 'data/callability/tracts/read_counts_tabix/{cross}.read_counts.tsv.gz'
    output:
        'data/rate/denominators/{cross}.tsv'
    params:
        window_size = '2000'
    shell:
        'time python analysis/rate/windowed_effective_sequence.py '
        ' --bam {input.bam} --snp_read_counts {input.tracts} '
        '--window_size {params.window_size} --out {output}'

rule windowed_snp_counts:
    input:
        vcf = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/rate/snp_density/{cross}.snps.2kb.tsv'
    params:
        window_size = '2000'
    script:
        'snp_density.py'

rule annotation_denominator_calc:
    input:
        bam = 'data/alignments/bam_prepped/sorted/{cross}.sorted.bam',
        gff_db = 'data/references/CC4532.v1_1.genes.primaryTrs.db',
        tracts = 'data/callability/tracts/read_counts_tabix/{cross}.read_counts.tsv.gz'
    output:
        'data/callability/callables/annotations/{cross}.tsv'
    shell:
        'time python analysis/rate/annotation_effective_sequence_simplified.py '
        '--bam {input.bam} --gff {input.gff_db} '
        '--snp_read_counts {input.tracts} --out {output}'
