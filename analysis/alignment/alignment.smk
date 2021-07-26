"""
alignment snakefile

usage:
snakemake -p \
--snakefile analysis/alignment/alignment.smk \
[rule]

"""

import os
from glob import glob

"""
SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/fastq/*x*.fastq.gz')]
PREFIXES = list(set([f.split('_')[0] for f in SAMPLES]))
TRIM_PREFIXES = [f'{p}_trim_{i}' for p in PREFIXES for i in [1, 2]]
"""
JAVA_EXEC = '/usr/bin/java'
PICARD = 'bin/picard.jar'

with open('data/alignments/samples.txt', 'r') as f:
    SAMPLES = [sample.rstrip() for sample in f]

rule all:
    input:
        expand("data/alignments/bam_temp/{sample}.fixMate.bam", sample=SAMPLES)

rule bwa_aln:
    input:
        ref = "data/references/CC4532.w_organelles_MTplus.fa",
        fwd = "data/alignments/fastq_trim/{sample}_trim_1.fq.gz",
        rev = "data/alignments/fastq_trim/{sample}_trim_2.fq.gz",
    output:
        temp("data/alignments/bam_temp/{sample}.bam")
    threads:
        16
    shell:
        "time bwa mem -t {threads} {input.ref} {input.fwd} {input.rev} | samtools view -Sb - > {output}"
        
rule bam_sort:
    input:
        bam = "data/alignments/bam_temp/{sample}.bam"
    output:
        "data/alignments/bam_temp/{sample}.sorted.bam"
    threads:
        4
    shell:
        "samtools sort -@{threads} -T {input.bam}.sorting.tmp "
        "-O bam -o {output} {input.bam}"
        
rule bam_fix_mate:
    input:
        bam = "data/alignments/bam_temp/{sample}.sorted.bam"
    output:
        "data/alignments/bam_temp/{sample}.fixMate.bam"
    shell:
        "{JAVA_EXEC} -jar {PICARD} FixMateInformation "
        "I={input.bam} O={output} "
        "VALIDATION_STRINGENCY=LENIENT"

