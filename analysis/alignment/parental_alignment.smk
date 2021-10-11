"""
parental alignment snakefile - generates parental bams 
for downstream variant calling (see `analysis/genotyping/variant_calling.smk`)

requires:
bwa, samtools, picard

usage:
snakemake -pr -s analysis/alignment/parental_alignment.smk 
"""

import os
from glob import glob

# --- globals

configfile: 'analysis/alignment/config.yml'
JAVA_EXEC = config['JAVA_EXEC']
PICARD = config['PICARD']

SAMPLES = [os.path.basename(f.rstrip('fastq.gz'))
    for f in glob('data/alignments/parental_fastq/*.fastq.gz')]
PREFIXES = list(set([re.search('(^[A-Z]{2}[0-9]{3,}[_590]*)_[12]', f).group(1)
    for f in SAMPLES]))

# handling CC2932 - there are two separate files
for prefix in ['CC2932_50', 'CC2932_90']:
    PREFIXES.remove(prefix)
PREFIXES.append('CC2932')

# --- rules

rule all:
    """
    the final outfiles are set as the index files, not the bams themselves,
    since indexing is the final step!
    """
    input:
        expand("data/alignments/parental_bam/{sample}.bam.bai", sample=PREFIXES)

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
        "time bwa mem -t {threads} {input.ref} {input.fwd} {input.rev} > {output}"

rule bam_convert:
    input:
        sam = "data/alignments/parental_bam_temp/{sample}.sam"
    output:
        temp("data/alignments/parental_bam_temp/{sample}.bam")
    shell:
        "samtools view -Sb {input.sam} > {output}"

rule combine_2932:
    """
    the 2932 reads were generated using teh a hybrid 2x50 and 2x90 approach.
    this rule combined the aligned reads into a single CC2932 bam
    """
    input:
        file_50 = "data/alignments/parental_bam_temp/CC2932_50.bam",
        file_90 = "data/alignments/parental_bam_temp/CC2932_90.bam"
    output:
        temp("data/alignments/parental_bam_temp/CC2932.bam")
    run:
        shell("echo {input.file_50} > temp_cat_files.txt")
        shell("echo {input.file_90} >> temp_cat_files.txt")
        shell("samtools cat -b temp_cat_files.txt -o {output}")
        shell("rm -v temp_cat_files.txt")
        
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
        "data/alignments/parental_bam/{sample}.bam"
    output:
        "data/alignments/parental_bam/{sample}.bam.bai"
    shell:
        "samtools index {input} {output}"
