"""
parental alignment snakefile

usage:
snakemake -p \
--snakefile analysis/alignment/parental_alignment.smk \
[rule]

"""

import os
from glob import glob

# globals
configfile: 'analysis/alignment/config.yml'
JAVA_EXEC = config['JAVA_EXEC']
PICARD = config['PICARD']

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/parental_fastq/*.fastq.gz')]
PREFIXES = list(set([re.search('(^[A-Z]{2}[0-9]{3,}[_590]*)_[12]', f).group(1)
    for f in SAMPLES]))

# functions
rule all:
    input:
        expand("data/alignments/parental_bam/{sample}.bam", sample=PREFIXES)

rule bwa_aln:
    input:
        ref = "data/references/CC4532.w_organelles_MTplus.fa",
        fwd = "data/alignments/parental_fastq_trim/{sample}_trim_1.fq.gz",
        rev = "data/alignments/parental_fastq_trim/{sample}_trim_2.fq.gz",
    output:
        temp("data/alignments/parental_bam_temp/{sample}.sam")
    threads:
        20
    shell:
        # "time bwa mem -t {threads} {input.ref} {input.fwd} {input.rev} | samtools view -Sb - > {output}"
        "time bwa mem -t {threads} {input.ref} {input.fwd} {input.rev} > {output}"

rule bam_convert:
    input:
        sam = "data/alignments/parental_bam_temp/{sample}.sam"
    output:
        temp("data/alignments/parental_bam_temp/{sample}.bam")
    shell:
        "samtools view -Sb {input.sam} > {output}"
        
rule bam_sort:
    input:
        bam = "data/alignments/parental_bam_temp/{sample}.bam"
    output:
        temp("data/alignments/parental_bam_temp/{sample}.sorted.bam")
    threads:
        4
    shell:
        "samtools sort -@{threads} -T {input.bam}.sorting.tmp "
        "-O bam -o {output} {input.bam}"
        
rule bam_fix_mate:
    input:
        bam = "data/alignments/parental_bam_temp/{sample}.sorted.bam"
    output:
        temp("data/alignments/parental_bam_temp/{sample}.fixMate.bam")
    shell:
        "{JAVA_EXEC} -jar {PICARD} FixMateInformation "
        "I={input.bam} O={output} "
        "VALIDATION_STRINGENCY=LENIENT"

rule bam_read_groups:
    input:
        bam = "data/alignments/parental_bam_temp/{sample}.fixMate.bam"
    output:
        temp("data/alignments/parental_bam_temp/{sample}.RG.bam")
    run:
        sample_name = os.path.basename(input.bam).rstrip('fixMate.RG.bam')
        shell("{JAVA_EXEC} -jar {PICARD} AddOrReplaceReadGroups "
              "I={input.bam} O={output} "
              "RGID={sample_name} RGLB=lib1 RGPL=illumina "
              "RGPU=unit1 RGSM={sample_name} "
              "VALIDATION_STRINGENCY=LENIENT")

rule bam_mark_duplicates:
    input:
        bam = "data/alignments/parental_bam_temp/{sample}.RG.bam"
    output:
        bam_out = "data/alignments/parental_bam/{sample}.bam",
        metrics = "data/alignments/parental_bam_temp/{sample}_dup_metrics.txt"
    shell:
        "{JAVA_EXEC} -jar {PICARD} MarkDuplicates "
        "I={input.bam} O={output.bam_out} "
        "REMOVE_DUPLICATES=true METRICS_FILE={output.metrics}"

rule bam_idx:
    input:
        expand("data/alignments/parental_bam/{sample}.bam", sample=PREFIXES)
    log: 
        "data/alignments/parental_bam/idx.log"
    run:
        for fname in input:
            shell(f"samtools index {fname}")
