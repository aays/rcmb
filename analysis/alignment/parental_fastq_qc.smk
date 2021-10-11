"""
parental read qc snakefile - for QC'ing and trimming parental reads

requires:
fastqc, trimmomatic

usage:
snakemake -pr -s analysis/alignment/parental_fastq_qc.smk [rul]

no `rule all` - needs to be run with explicit rules at the command line
"""

import re
import os
from glob import glob

# --- globals

configfile: 'analysis/alignment/config.yml'

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/parental_fastq/*.fastq.gz')]
PREFIXES = list(set([re.search('(^[A-Z]{2}[0-9]{3,}[_590]*)_[12]', f).group(1)
    for f in SAMPLES]))
TRIM_PREFIXES = [f'{p}_trim_{i}' for p in PREFIXES for i in [1,2]]

# --- rules

rule fastqc:
    input:
        expand("data/alignments/parental_fastq/{sample}.fastq.gz", sample=SAMPLES)
    log:
        "data/alignments/parental_fastq/fastqc/fastqc.log"
    shell:
        "fastqc --kmers 7 --outdir data/alignments/parental_fastq/fastqc/raw {input}"

rule unzip_fastqc:
    input:
        expand("data/alignments/parental_fastq/fastqc/raw/{sample}_fastqc.zip", sample=SAMPLES)
    log:
        "data/alignments/parental_fastq/fastqc/fastqc.log"
    run:
        for fname in input:
            shell(f"unzip {fname} -d data/alignments/parental_fastq/fastqc/")

rule trim_reads:
    input:
        expand("data/alignments/parental_fastq/{sample}.fastq.gz", sample=SAMPLES)
    log:
        "data/alignments/parental_fastq_trim/trim.log"
    threads:
        16
    run:
        for prefix in PREFIXES:
            shell('trimmomatic PE -threads {threads} -phred33 '
                  'data/alignments/parental_fastq/{prefix}_1.fastq.gz '
                  'data/alignments/parental_fastq/{prefix}_2.fastq.gz '
                  'data/alignments/parental_fastq_trim/{prefix}_trim_1.fq.gz '
                  'data/alignments/parental_fastq_trim/{prefix}_trim_unpaired_1.fq.gz '
                  'data/alignments/parental_fastq_trim/{prefix}_trim_2.fq.gz '
                  'data/alignments/parental_fastq_trim/{prefix}_trim_unpaired_2.fq.gz '
                  'SLIDINGWINDOW:5:30 LEADING:5 TRAILING:5 2>> {log}')

rule fastqc_trim:
    input:
        expand("data/alignments/parental_fastq_trim/{trim_file}.fq.gz",
            trim_file=TRIM_PREFIXES)
    log:
        "data/alignments/parental_fastq_trim/fastqc/fastqc.log"
    shell:
        "fastqc --kmers 7 --outdir data/alignments/parental_fastq_trim/fastqc/raw {input}"

rule unzip_fastqc_trim:
    input:
        expand("data/alignments/parental_fastq_trim/fastqc/raw/{trim_file}_fastqc.zip",
            trim_file=sorted(TRIM_PREFIXES))
    log:
        "data/alignments/parental_fastq_trim/fastqc/fastqc.log"
    run:
        for fname in input:
            shell(f"unzip {fname} -d data/alignments/parental_fastq_trim/fastqc/")
