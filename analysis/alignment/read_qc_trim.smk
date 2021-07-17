"""
read qc snakefile

usage:
snakemake -p \
--snakefile analysis/alignment/read_qc_trim.smk \
[rule]

"""

import os
from glob import glob

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/fastq/*x*.fastq.gz')]

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
