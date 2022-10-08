"""
ldhelmet workflow - generates per-sample fastas from previously
generated GVCFs, converts them to per-chr fastas, and then runs
LDhelmet on them

requires:
bcftools, gatk4, ldhelmet

usage:
snakemake -pr -s analysis/ldhelmet/ldhelmet.smk --cores [cores]
"""

import os
from glob import glob

# --- globals

SAMPLES = list(set([
    os.path.basename(fname.rstrip('.bam'))
    for fname in glob('data/alignments/parental_bam/*')
    if fname.endswith('.bam')]
))

CHROM_LENGTHS = {
    'chromosome_01': 8225636,
    'chromosome_02': 8655884,
    'chromosome_03': 9286894,
    'chromosome_04': 4130073,
    'chromosome_05': 3682160,
    'chromosome_06': 8913359,
    'chromosome_07': 6492107,
    'chromosome_08': 4526983,
    'chromosome_09': 6807148,
    'chromosome_10': 6800247,
    'chromosome_11': 4479522,
    'chromosome_12': 9952739,
    'chromosome_13': 5281438,
    'chromosome_14': 4217303,
    'chromosome_15': 5870643,
    'chromosome_16': 8042475,
    'chromosome_17': 6954842}


# --- paths

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        print(f'[rcmb] {dir_path} exists')

mkdir('data/ldhelmet/fasta')
mkdir('data/ldhelmet/fasta_chrom')

# --- rules

rule all:
    input:
        expand('data/ldhelmet/fasta/{sample}.{chrom}.fa', 
            sample=SAMPLES, chrom=list(CHROM_LENGTHS.keys())),
        expand('data/ldhelmet/fasta_chrom/{pop}.chromosome_0{n}.fa', 
            n=list(range(1, 10)), pop=['NA1', 'NA2', 'ALL']),
        expand('data/ldhelmet/fasta_chrom/{pop}.chromosome_{n}.fa', 
            n=list(range(10, 18)), pop=['NA1', 'NA2', 'ALL']),
        expand('data/ldhelmet/{pop}/chromosome_0{n}.txt', 
            n=list(range(1, 10)), pop=['NA1', 'NA2']),
        expand('data/ldhelmet/{pop}/chromosome_{n}.txt', 
            n=list(range(10, 18)), pop=['NA1', 'NA2'])

rule create_indiv_fastas:
    input:
        ref = 'data/references/CC4532.w_organelles_MTplus.fa'
    output:
        expand(
            'data/ldhelmet/fasta/{sample}.{chrom}.fa',
            sample=SAMPLES, chrom=list(CHROM_LENGTHS.keys()))
    params:
        min_GQ = '30'
    run:
        for n in range(1, 18):
            chrom_name = f'chromosome_0{n}' if n < 10 else f'chromosome_{n}'
            chrom_length = CHROM_LENGTHS[chrom_name]
            for sample in SAMPLES:
                gvcf_path = f'data/genotyping/gvcf_sample/{sample}.g.vcf.gz'
                shell(
                    'time python bin/vcf2fasta.py -v {gvcf_path} '
                    '-r {input.ref} -s {sample} -i {chrom_name}:1-{chrom_length} '
                    '--min_GQ {params.min_GQ} > data/ldhelmet/fasta/{sample}.{chrom_name}.fa'
                )
        

rule create_per_chr_fastas:
    input:
        expand('data/ldhelmet/fasta/{sample}.{chrom}.fa', 
            sample=SAMPLES, chrom=list(CHROM_LENGTHS.keys()))
    output:
        expand('data/ldhelmet/fasta_chrom/ALL.chromosome_0{n}.fa', n=list(range(1, 10))),
        expand('data/ldhelmet/fasta_chrom/NA1.chromosome_0{n}.fa', n=list(range(1, 10))),
        expand('data/ldhelmet/fasta_chrom/NA2.chromosome_0{n}.fa', n=list(range(1, 10))),
        expand('data/ldhelmet/fasta_chrom/ALL.chromosome_{n}.fa', n=list(range(10, 18))),
        expand('data/ldhelmet/fasta_chrom/NA1.chromosome_{n}.fa', n=list(range(10, 18))),
        expand('data/ldhelmet/fasta_chrom/NA2.chromosome_{n}.fa', n=list(range(10, 18)))
    script:
        'convert_fastas.py'

rule ldhelmet_find_confs:
    input:
        'data/ldhelmet/fasta_chrom/{pop}.{chrom}.fa'
    output:
        'data/ldhelmet/{pop}/{chrom}.conf'
    params:
        window_size = 50
    threads:
        10
    shell:
        'time ./bin/ldhelmet find_confs --num_threads {threads} '
        '--window_size {params.window_size} '
        '--output_file {output} {input}'

rule ldhelmet_table_gen:
    input:
        'data/ldhelmet/{pop}/{chrom}.conf'
    output:
        'data/ldhelmet/{pop}/{chrom}.lk'
    log:
        'data/ldhelmet/{pop}/table_gen.log'
    threads:
        10
    params:
        theta = '0.03',
        rhos = '0.0 0.1 10.0 1.0 100.0',
        log_2 = 'data/ldhelmet/{pop}/table_gen_2.log'
    shell:
        'time ./bin/ldhelmet table_gen --num_threads {threads} '
        '--conf_file {input} --theta {params.theta} '
        '--rhos {params.rhos} --output_file {output} > {log} 2> {params.log_2}'

rule ldhelmet_pade:
    input:
        temp('data/ldhelmet/{pop}/{chrom}.conf')
    output:
        'data/ldhelmet/{pop}/{chrom}.pade'
    params:
        theta = '0.03'
    shell:
        'time ./bin/ldhelmet pade --num_threads {threads} '
        '--conf_file {input} --theta {params.theta} '
        '--output_file {output}'

rule ldhelmet_rjmcmc:
    input:
        seq_file = 'data/ldhelmet/fasta_chrom/{pop}.{chrom}.fa',
        lk_file = temp('data/ldhelmet/{pop}/{chrom}.lk'),
        pade_file = temp('data/ldhelmet/{pop}/{chrom}.pade'),
        mut_mat_file = 'data/ldhelmet/mut_mat'
    output:
        'data/ldhelmet/{pop}/{chrom}.post'
    threads:
        30
    params:
        window_size = 50,
        num_iter = '1000000',
        burn_in = '100000',
        block_penalty = '100'
    shell:
        'time ./bin/ldhelmet rjmcmc '
        '--num_threads {threads} --seq_file {input.seq_file} '
        '--lk_file {input.lk_file} --pade_file {input.pade_file} '
        '--num_iter {params.num_iter} --burn_in {params.burn_in} '
        '--block_penalty {params.block_penalty} '
        '--mut_mat_file {input.mut_mat_file} '
        '--output_file {output}'

rule ldhelmet_post_to_text:
    input:
        'data/ldhelmet/{pop}/{chrom}.post'
    output:
        'data/ldhelmet/{pop}/{chrom}.txt'
    shell:
        'time ./bin/ldhelmet post_to_text '
        '--mean --perc 0.025 --perc 0.50 --perc 0.975 '
        '--output_file {output} {input}'


    
        
        


