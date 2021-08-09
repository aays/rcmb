"""
read qc snakefile

usage:
snakemake -p \
--snakefile analysis/alignment/read_qc_trim.smk \
[rule]

"""

import os
from glob import glob

# globals
configfile: 'config.yml'

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/fastq/*x*.fastq.gz')]
PREFIXES = list(set([f.split('_')[0] for f in SAMPLES]))
TRIM_PREFIXES = [f'{p}_trim_{i}' for p in PREFIXES for i in [1, 2]]
JAVA_EXEC = config['JAVA_EXEC']

# functions
rule fastqc:
    input:
        expand("data/alignments/fastq/{sample}.fastq.gz", sample=SAMPLES)
    log:
        "data/alignments/fastq/fastqc/fastqc.log"
    shell:
        "fastqc --kmers 7 --outdir data/alignments/fastq/fastqc/raw {input}"

rule unzip_fastqc:
    input:
        expand("data/alignments/fastq/fastqc/raw/{sample}_fastqc.zip", sample=SAMPLES)
    log:
        "data/alignments/fastq/fastqc/fastqc.log"
    run:
        for fname in input:
            shell(f"unzip {fname} -d data/alignments/fastq/fastqc/")

rule trim_reads:
    input:
        expand("data/alignments/fastq/{sample}.fastq.gz", sample=SAMPLES)
    log:
        "data/alignments/fastq_trim/trim.log"
    threads:
        16
    run:
        for prefix in PREFIXES:
            shell('trimmomatic PE -threads {threads} -phred33 '
                  'data/alignments/fastq/{prefix}_R1.fastq.gz '
                  'data/alignments/fastq/{prefix}_R2.fastq.gz '
                  'data/alignments/fastq_trim/{prefix}_trim_1.fq.gz '
                  'data/alignments/fastq_trim/{prefix}_trim_unpaired_1.fq.gz '
                  'data/alignments/fastq_trim/{prefix}_trim_2.fq.gz '
                  'data/alignments/fastq_trim/{prefix}_trim_unpaired_2.fq.gz '
                  'ILLUMINACLIP:bin/NEBNext_dual.fasta:2:30:10 '
                  'SLIDINGWINDOW:4:20')

rule fastqc_trim:
    input:
        expand("data/alignments/fastq_trim/{trim_file}.fq.gz", 
            trim_file=TRIM_PREFIXES)
    log:
        "data/alignments/fastq_trim/fastqc/fastqc.log"
    shell:
        "fastqc --kmers 7 --outdir data/alignments/fastq_trim/fastqc/raw {input}"

rule unzip_fastqc_trim:
    input:
        expand("data/alignments/fastq_trim/fastqc/raw/{trim_file}_fastqc.zip", 
            trim_file=sorted(TRIM_PREFIXES))
    log:
        "data/alignments/fastq_trim/fastqc/fastqc.log"
    run:
        for fname in input:
            shell(f"unzip {fname} -d data/alignments/fastq_trim/fastqc/")


