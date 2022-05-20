"""
variant calling workflow - generates parental pair VCFs for
use in recombination event detection

redoing with freebayes

requires:
freebayes, readcomb

usage:
snakemake -pr -s analysis/genotyping/variant_calling.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

SAMPLES = list(set([
    os.path.basename(fname.rstrip('.bam'))
    for fname in glob('data/alignments/parental_bam/*')
    if fname.endswith('.bam')]
))

with open('data/genotyping/samples.txt', 'r') as f:
    CROSSES = [line.split(' ')[0] for line in f]
    CROSS_DICT = {}
    for line in CROSSES:
        mt_plus, mt_minus = line.split('x')
        if mt_plus not in CROSS_DICT:
            CROSS_DICT[mt_plus] = [mt_minus]
        else:
            CROSS_DICT[mt_plus].append(mt_minus)

CHROMS = [f'chromosome_0{i}' for i in range(1, 10)]
CHROMS.extend([f'chromosome_{i}' for i in range(10, 18)])
CHROMS.extend(['cpDNA', 'mtDNA'])

CONCAT_DICT = {}
for cross in CROSSES:
    CONCAT_DICT[cross] = []
    for chrom in CHROMS:
        CONCAT_DICT[cross].append(f'data/genotyping/vcf_chrom/{cross}.{chrom}.vcf.gz')

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/genotyping/vcf_chrom')
mkdir('data/genotyping/vcf_freebayes')
mkdir('data/genotyping/vcf_filtered')

# --- rules

rule all:
    input:
        # expand('data/genotyping/vcf_chrom/{cross}.{chrom}.vcf', cross=CROSSES, chrom=CHROMS),
        expand('data/genotyping/vcf_freebayes/{cross}.vcf.gz', cross=CROSSES),
        expand('data/genotyping/vcf_filtered/{cross}.vcf.gz.tbi', cross=CROSSES)

rule freebayes_chrom:
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa'
    output:
        expand('data/genotyping/vcf_chrom/{cross}.{chrom}.vcf', cross=CROSSES, chrom=CHROMS)
    threads:
        6
    params:
        ploidy = '2'
    run:
        for cross in CROSSES:
            mt_plus, mt_minus = [
                f'data/alignments/parental_bam/CC{sample}.bam' for sample in cross.split('x')]
            mt_plus = mt_plus.replace('CCGB', 'GB') # awkward - handles GB119

            # if cross has already been called - skip over
            if len(glob(f'data/genotyping/vcf_chrom/{cross}*')) > 0:
                print(f'[rcmb] skipping {cross} - already completed')
                continue
            else:
                shell(
                    "parallel -j {threads} "
                    "'freebayes -f {input.ref} --theta 0.02 --ploidy {params.ploidy} "
                    "-r chromosome_{{}} --genotype-qualities "
                    "{mt_plus} {mt_minus} > data/genotyping/vcf_chrom/{cross}.chromosome_{{}}.vcf'"
                    " ::: {{01..09}} {{10..17}}"
                )
                # organelles
                shell(
                    "parallel -j {threads} "
                    "'freebayes -f {input.ref} --theta 0.02 --ploidy {params.ploidy} "
                    "-r {{}} --genotype-qualities "
                    "{mt_plus} {mt_minus} > data/genotyping/vcf_chrom/{cross}.{{}}.vcf'"
                    " ::: mtDNA cpDNA"
                )
                shell('sleep 1')
                shell('cp -v data/genotyping/vcf_chrom/{cross}* data/genotyping/vcf_chrom_backup')

rule bgzip_tabix_chroms:
    input:
        vcf = 'data/genotyping/vcf_chrom/{cross_chrom}.vcf'
    output:
        'data/genotyping/vcf_chrom/{cross_chrom}.vcf.gz'
    shell:
        "bgzip {input.vcf}; tabix {input.vcf}.gz"
        

rule vcf_concat:
    input:
        vcf_in = temp(lambda wildcards: CONCAT_DICT[wildcards.cross])
    output:
        vcf_out = 'data/genotyping/vcf_freebayes/{cross}.vcf'
    threads:
        6
    shell:
        "bcftools concat {input.vcf_in} > {output.vcf_out}"


rule bgzip_tabix_combined:
    input:
        vcf = 'data/genotyping/vcf_freebayes/{cross}.vcf'
    output:
        'data/genotyping/vcf_freebayes/{cross}.vcf.gz'
    shell:
        "bgzip {input.vcf}; tabix {input.vcf}.gz"
    

rule readcomb_vcfprep:
    input:
        vcf_in = 'data/genotyping/vcf_freebayes/{cross}.vcf.gz',
    output:
        'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    params:
        min_qual = '30',
        purity_filter = '1'
    shell:
        'time readcomb-vcfprep --vcf {input.vcf_in} --no_hets --snps_only '
        '--min_GQ {params.min_qual} --purity_filter {params.purity_filter} --out {output}'

rule vcf_tabix:
    input:
        vcf_in = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/genotyping/vcf_filtered/{cross}.vcf.gz.tbi'
    shell:
        'tabix -p vcf {input.vcf_in}'

