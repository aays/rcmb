"""
phase change detection workflow

requires:
samtools, readcomb

usage:
snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

with open('data/genotyping/samples.txt', 'r') as f:
    CROSSES = [line.split(' ')[0] for line in f]
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

mkdir('data/phase_changes/sam')
mkdir('data/phase_changes/bam')
mkdir('data/phase_changes/fp_bed')

# --- rules

rule all:
    input:
        expand('data/phase_changes/sam/{cross}.init.sam', cross=CROSSES),
        expand('data/phase_changes/bam/{cross}.filtered.bam', cross=CROSSES),
        expand('data/phase_changes/bam/{cross}.filtered.bam.bai', cross=CROSSES),
        expand('data/phase_changes/event_summaries/{cross}.75.tsv', cross=CROSSES)

rule readcomb_filter:
    input:
        bam = 'data/alignments/bam_prepped/{cross}.sorted.bam',
        vcf = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/phase_changes/sam/{cross}.init.sam'
    threads:
        24
    params:
        log = 'data/phase_changes/phase_change_filter.log',
        min_qual = '0',
        min_mapq = '0'
    shell:
        'time readcomb-filter --bam {input.bam} --vcf {input.vcf} '
        '--processes {threads} --log {params.log} --quality {params.min_qual} '
        '--min_mapq {params.min_mapq} --out {output}'

rule readcomb_fp:
    input:
        sam = 'data/phase_changes/sam/{cross}.init.sam',
        false_plus = 'data/phase_changes/parental/{cross}.plus.filtered.sam',
        false_minus = 'data/phase_changes/parental/{cross}.minus.filtered.sam',
        vcf = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        bed_lookup = 'data/phase_changes/fp_bed/{cross}.midpoint.bed.gz',
        filtered_sam = 'data/phase_changes/sam/{cross}.filtered.sam'
    params:
        log = 'data/phase_changes/false_positives.log',
        method = 'midpoint',
        length = '0', # read length filter
        qual = '0' # base qual
    shell:
        'time readcomb-fp --fname {input.sam} '
        '--false_plus {input.false_plus} --false_minus {input.false_minus} '
        '--vcf {input.vcf} --method {params.method} --read_length_filter {params.length} '
        '--base_qual_filter {params.qual} --false_bed_out {output.bed_lookup} '
        '--log {params.log} --out {output.filtered_sam}'
    
rule bam_sort:
    input:
        sam = 'data/phase_changes/sam/{cross}.filtered.sam'
    output:
        'data/phase_changes/bam/{cross}.filtered.bam'
    shell:
        'samtools sort -O bam {input.sam} > {output}'

rule bam_idx:
    input:
        bam = 'data/phase_changes/bam/{cross}.filtered.bam'
    output:
        'data/phase_changes/bam/{cross}.filtered.bam.bai'
    shell:
        'samtools index {input.bam} {output}'

rule summarise_cross:
    input:
        bam = 'data/phase_changes/sam/{cross}.filtered.sam',
        vcf = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/phase_changes/event_summaries/{cross}.75.tsv'
    threads:
        12
    params:
        mapq = '0',
        mask_size = '75',
        base_qual = '0'
    shell:
        'python analysis/phase_changes/summarise_cross.py '
        '--bam {input.bam} --vcf {input.vcf} --mask_size {params.mask_size} '
        '--processes {threads} --mapq {params.mapq} --base_qual {params.base_qual} ' 
        '--out {output}'

