"""
cross alignment snakefile - generates cross bams for direct
recombination event detection

requires:
bwa, samtools, picard

usage:
snakemake -pr -s analysis/alignment/alignment.smk

"""

import os
from glob import glob

# --- globals

configfile: 'analysis/alignment/config.yml'
JAVA_EXEC = config['JAVA_EXEC']
PICARD = config['PICARD']

with open('data/alignments/samples.txt', 'r') as f:
    SAMPLES = [sample.rstrip() for sample in f]

# --- rules

rule all:
    input:
        expand("data/alignments/bam/{sample}.bam.bai", sample=SAMPLES)

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
        temp("data/alignments/bam_temp/{sample}.sorted.bam")
    threads:
        4
    shell:
        "samtools sort -@{threads} -T {input.bam}.sorting.tmp "
        "-O bam -o {output} {input.bam}"
        
rule bam_fix_mate:
    input:
        bam = "data/alignments/bam_temp/{sample}.sorted.bam"
    output:
        temp("data/alignments/bam_temp/{sample}.fixMate.bam")
    shell:
        "{JAVA_EXEC} -jar {PICARD} FixMateInformation "
        "I={input.bam} O={output} "
        "VALIDATION_STRINGENCY=LENIENT"

rule bam_read_groups:
    input:
        bam = "data/alignments/bam_temp/{sample}.fixMate.bam"
    output:
        "data/alignments/bam_temp/{sample}.RG.bam"
    run:
        sample_name = os.path.basename(input.bam).rstrip('fixMate.RG.bam')
        shell("{JAVA_EXEC} -jar {PICARD} AddOrReplaceReadGroups "
              "I={input.bam} O={output} "
              "RGID={sample_name} RGLB=lib1 RGPL=illumina "
              "RGPU=unit1 RGSM={sample_name} "
              "VALIDATION_STRINGENCY=LENIENT")

rule bam_mark_duplicates:
    input:
        bam = "data/alignments/bam_temp/{sample}.RG.bam"
    output:
        bam_out = "data/alignments/bam/{sample}.bam",
        metrics = "data/alignments/bam_temp/{sample}_dup_metrics.txt"
    shell:
        "{JAVA_EXEC} -jar {PICARD} MarkDuplicates "
        "I={input.bam} O={output.bam_out} "
        "REMOVE_DUPLICATES=true METRICS_FILE={output.metrics}"

rule bam_idx:
    input:
        bam = "data/alignments/bam/{sample}.bam"
    output:
        "data/alignments/bam/{sample}.bam.bai"
    shell:
        "samtools index {input.bam} {output}"
