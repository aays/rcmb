"""
parental phase change detection - used to curate false positives that
can then be removed from the recombinant read set

requires:
samtools, readcomb

usage:
snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

with open('data/genotyping/samples.txt', 'r') as f:
    CROSSES = [line.split(' ')[0] for line in f]
    CROSS_DICT = {}
    PARENT_CROSSES = []
    STRAINS = {'mt_plus': [], 'mt_minus': []}
    for line in CROSSES:
        mt_plus, mt_minus = line.split('x')
        PARENT_CROSSES.append(f'{mt_plus}x{mt_minus}.plus')
        PARENT_CROSSES.append(f'{mt_plus}x{mt_minus}.minus')
        STRAINS['mt_plus'].append(mt_plus)
        STRAINS['mt_minus'].append(mt_minus)
        if mt_plus not in CROSS_DICT:
            CROSS_DICT[mt_plus] = [mt_minus]
        else:
            CROSS_DICT[mt_plus].append(mt_minus)


STRAINS['mt_plus'] = list(set(STRAINS['mt_plus']))
STRAINS['mt_minus'] = list(set(STRAINS['mt_minus']))
PARENTS = STRAINS['mt_plus'] + STRAINS['mt_minus']
PARENTS = ['CC' + parent for parent in PARENTS if not parent.startswith('GB')]

MINUS_CROSS_DICT = {}
for mt_minus in STRAINS['mt_minus']:
    MINUS_CROSS_DICT[mt_minus] = [mt_plus for mt_plus in CROSS_DICT if mt_minus in CROSS_DICT[mt_plus]]

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/alignments/parental_bam_filtered')
mkdir('data/phase_changes/parental')

# --- rules

rule all:
    input:
        expand('data/alignments/parental_bam_filtered/{parent}.sorted.bam', parent=PARENTS)

rule readcomb_bamprep:
    input:
        bam = "data/alignments/parental_bam/{parent}.bam"
    output:
        "data/alignments/parental_bam_filtered/{parent}.sorted.bam"
    threads:
        16
    shell:
        "time readcomb-bamprep --bam {input.bam} "
        "--threads {threads} --outdir data/alignments/parental_bam_filtered/"

rule readcomb_filter:
    input:
        expand('data/alignments/parental_bam_filtered/{parent}.sorted.bam', parent=PARENTS)
    threads:
        16
    params:
        log = 'data/phase_changes/parental_phase_changes.log',
        min_qual = '30'
    run:
        for mt_plus in CROSS_DICT: 
            cross_strains = CROSS_DICT[mt_plus]
            for mt_minus in cross_strains:
                cross = f'{mt_plus}x{mt_minus}'
                if not mt_plus.startswith('GB') and not mt_plus.startswith('CC'):
                    mt_plus = 'CC' + mt_plus
                if 'CC' in cross:
                    cross = cross.replace('CC', '') # weird CC name handling

                # only run rule if file doesn't already exist
                if os.path.exists(f'data/phase_changes/parental/{cross}.plus.filtered.sam'):
                    print(f'[rcmb] {cross}.plus.filtered.sam exists. skipping...')
                    continue

                shell(
                    'time readcomb-filter --bam data/alignments/parental_bam_filtered/{mt_plus}.sorted.bam '
                    '--vcf data/genotyping/vcf_filtered/{cross}.vcf.gz ' 
                    '--processes {threads} --log {params.log} '
                    '--quality {params.min_qual} '
                    '--out data/phase_changes/parental/{cross}.plus.filtered.sam')
        for mt_minus in MINUS_CROSS_DICT:
            cross_strains = MINUS_CROSS_DICT[mt_minus]
            for mt_plus in cross_strains:
                cross = f'{mt_plus}x{mt_minus}'
                if not mt_minus.startswith('GB') and not mt_minus.startswith('CC'):
                    minus = 'CC' + minus
                if 'CC' in cross:
                    cross = cross.replace('CC', '') # weird CC name handling

                # only run rule if file doesn't already exist
                if os.path.exists(f'data/phase_changes/parental/{cross}.minus.filtered.sam'):
                    print(f'[rcmb] {cross}.minus.filtered.sam exists. skipping...')
                    continue

                shell(
                    'time readcomb-filter --bam data/alignments/parental_bam_filtered/{mt_minus}.sorted.bam '
                    '--vcf data/genotyping/vcf_filtered/{cross}.vcf.gz '
                    '--processes {threads} --log {params.log} '
                    '--quality {params.min_qual} '
                    '--out data/phase_changes/parental/{cross}.minus.filtered.sam')

