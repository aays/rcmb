"""
callability workflow

requires:
cyvcf2

usage:
snakemake -pr -s analysis/callability/callability.smk --cores [cores]
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

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/parental_fastq/*.fastq.gz')]
PARENTS = list(set([re.search('[A-Z]{2}[0-9]{3,4}', f).group(0)
    for f in SAMPLES]))

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/callability/tracts/')

# --- rules

rule all:
    input:
        expand('data/callability/tracts/{cross}.tracts.tsv', cross=CROSSES),
        'data/callability/tract_distribution.tsv',
        expand('data/callability/coverage/{parent}.cov.temp.tsv', parent=PARENTS),
        expand('data/callability/tracts/{cross}.read_counts.tsv', cross=CROSSES),
        'data/callability/callables/chrom_callables.tsv',
        expand('data/callability/eff_bp/{cross}.chrom.tsv', cross=CROSSES)


rule create_snp_tract_bed:
    input:
        'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/callability/tracts/{cross}.tracts.tsv'
    script:
        'create_snp_tract_bed.py'

rule summarise_tracts:
    input:
        expand('data/callability/tracts/{cross}.tracts.tsv', cross=CROSSES)
    output:
        'data/callability/tract_distribution.tsv'
    script:
        'summarise_tracts.py'

rule windowed_coverage:
    input:
        'data/alignments/parental_bam/{parent}.bam'
    output:
        'data/callability/coverage/{parent}.cov.temp.tsv'
    params:
        window_size = 2000
    script:
        'windowed_coverage.py'

rule get_tract_counts:
    input:
        tracts = 'data/callability/tracts/{cross}.tracts.tsv',
        bam = 'data/alignments/bam_prepped/sorted/{cross}.sorted.bam'
    output:
        'data/callability/tracts/{cross}.read_counts.tsv'
    script:
        'get_tract_read_counts.py'

rule summarise_chroms:
    input:
        expand('data/callability/tracts/{cross}.read_counts.tsv', cross=CROSSES)
    output:
        'data/callability/callables/chrom_callables.tsv'
    script:
        'summarise_chroms.py'

rule summarise_chrom_eff_bp:
    input:
        'data/callability/tracts/{cross}.read_counts.tsv'
    output:
        'data/callability/eff_bp/{cross}.chrom.tsv'
    script:
        'summarise_chrom_eff_bp.py'

    
