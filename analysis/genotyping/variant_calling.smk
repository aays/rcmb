"""
variant calling workflow - generates parental pair VCFs for
use in recombination event detection

follows the GATK germline short variant discovery
pipeline (per sample calls -> consolidation -> joint genotyping)
before filtering VCFs with readcomb-vcfprep

requires:
gatk4, readcomb

usage:
snakemake -pr -s analysis/genotyping/variant_calling.smk --cores [cores]
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

# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/genotyping/gvcf_sample')
mkdir('data/genotyping/gvcf_combined')
mkdir('data/genotyping/vcf')

# --- rules

rule all:
    input:
        expand('data/genotyping/vcf_filtered/{cross}.vcf.gz.tbi', cross=CROSSES)

rule call_indiv_gvcf:
    """
    recommend providing more cores at the command line (e.g. with --cores)
    to speed this up based if there's enough available memory - this will
    allow multiple GVCFs to be generated at once, even though each individual
    GVCF won't be generated any faster
    """
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa',
        bam = 'data/alignments/parental_bam/{sample}.bam'
    output:
        'data/genotyping/gvcf_sample/{sample}.g.vcf.gz'
    threads:
        2
    shell:
        'time gatk HaplotypeCaller -R {input.ref} '
        '-I {input.bam} -ERC GVCF -ploidy 2 '
        '--heterozygosity 0.02 --indel-heterozygosity 0.002 '
        '--native-pair-hmm-threads {threads} -O {output}'

rule combine_parents:
    """
    creates parental pair GVCFs - formats path while being
    mindful of CC/GB prefixes

    GB was included in the sample names in the cross bams but CC was not,
    but at the same time CC *was* a part of the GVCF names
    """
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa'
    output:
        expand('data/genotyping/gvcf_combined/{cross}.g.vcf.gz', cross=CROSSES)
    run:
        for mt_plus in CROSS_DICT:
            for mt_minus in CROSS_DICT[mt_plus]:

                if not mt_plus.startswith('GB'):
                    mt_plus_path = f'data/genotyping/gvcf_sample/CC{mt_plus}.g.vcf.gz'
                else:
                    mt_plus_path = f'data/genotyping/gvcf_sample/{mt_plus}.g.vcf.gz'

                if not mt_minus.startswith('GB'):
                    mt_minus_path = f'data/genotyping/gvcf_sample/CC{mt_minus}.g.vcf.gz'
                else:
                    mt_minus_path = f'data/genotyping/gvcf_sample/{mt_minus}.g.vcf.gz'

                cross_name = f'{mt_plus}x{mt_minus}'

                if os.path.exists(f'data/genotyping/gvcf_combined/{cross_name}.g.vcf.gz'):
                    print(f'[rcmb] {cross_name}.g.vcf.gz exists. skipping...')
                    continue
                else:
                    print(os.path.exists(f'data/genotyping/gvcf_combined/{cross_name}.g.vcf.gz'))
                    shell("gatk CombineGVCFs -R {input.ref} "
                      "--variant {mt_plus_path} --variant {mt_minus_path} "
                      "-O data/genotyping/gvcf_combined/{cross_name}.g.vcf.gz")

rule genotype_gvcfs:
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa',
        gvcf = 'data/genotyping/gvcf_combined/{cross}.g.vcf.gz'
    output:
        'data/genotyping/vcf/{cross}.vcf.gz'
    threads:
        2
    shell:
        'time gatk GenotypeGVCFs -R {input.ref} '
        '-V {input.gvcf} -O {output}'

rule readcomb_vcfprep:
    input:
        vcf_in = 'data/genotyping/vcf/{cross}.vcf.gz',
    output:
        'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    params:
        min_qual = '30',
        purity_filter = '1'
    threads:
        2
    shell:
        'time readcomb-vcfprep --vcf {input.vcf_in} --no_hets --snps_only '
        '--min_GQ {params.min_qual} --purity_filter {params.purity_filter} --out {output}'

rule vcf_tabix:
    input:
        vcf_in = 'data/genotyping/vcf_filtered/{cross}.vcf.gz'
    output:
        'data/genotyping/vcf_filtered/{cross}.vcf.gz.tbi'
    shell:
        'time tabix -p vcf {input.vcf_in}'

