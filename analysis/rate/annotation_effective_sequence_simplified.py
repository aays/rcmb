"""
windowed_effective_sequence.py - get effective sequence for a given
cross but split into windows
"""

import os
import math
import csv
import argparse
import pysam
import gffutils
from tqdm import tqdm
import readcomb.classification as rc

def args():
    parser = argparse.ArgumentParser(
        description='get effective sequence for GFF features', 
        usage='python windowed_effective_sequence.py [options]')

    parser.add_argument('-b', '--bam', required=True,
        type=str, help='BAM file containing all reads considered for phase changes')
    parser.add_argument('-g', '--gff', required=True,
        type=str, help='DB file containing GFF data')
    parser.add_argument('-s', '--snp_read_counts', required=True,
        type=str, help='SNP tract read count file')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.bam, args.gff, args.snp_read_counts, args.out

lengths = {
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
    'chromosome_17': 6954842
    }

def get_eff_bp(bam_reader, chrom, start, end):
    eff_bp = 0
    for p_column in bam_reader.pileup(chrom, start, end, truncate=True):
        for p_read in p_column.pileups:
            eff_bp += 1
    return eff_bp

def denominator_counts(bam, gff, snp_read_counts, out, lengths):

    db = gffutils.FeatureDB(gff)
    annotations = ['CDS', 'five_prime_UTR', 'three_prime_UTR']

    cross = os.path.basename(bam).rstrip('.sorted.bam')
    bam_reader = pysam.AlignmentFile(bam, 'rb')
    snp_reader = pysam.TabixFile(snp_read_counts)

    with open(out, 'w') as f:
        fieldnames = [
            'cross', 'chrom', 'CDS_count', 'CDS_eff_bp', 'utr5_count', 'utr5_eff_bp',
            'utr3_count', 'utr3_eff_bp', 'intronic_count', 'intronic_eff_bp',
            'intergenic_count', 'intergenic_eff_bp']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for chrom in lengths:

            # instantiate dict keeping feature counts and eff bp counts
            feature_dict = {
                'CDS': 0, 'five_prime_UTR': 0, 'three_prime_UTR': 0, 
                'inter_CDS_CDS': 0, 'inter_gene_gene': 0} # introns and intergenic
            eff_bp_dict = {
                'CDS': 0, 'five_prime_UTR': 0, 'three_prime_UTR': 0, 
                'inter_CDS_CDS': 0, 'inter_gene_gene': 0}

            # iterate through genes for all but intergenic
            region = f'{chrom}:1-{lengths[chrom]}'
            for i, gene in tqdm(enumerate(
                db.features_of_type('gene', limit=region)), desc=f'{chrom} gene features'):

                children = [f for f in db.children(gene, featuretype=annotations)]
                coding_features = [f for f in children if f.featuretype == 'CDS']
                introns = [f for f in db.interfeatures(coding_features)
                    if f.end > f.start] # avoids 'negative space' introns from overlapping features
                all_gene_features = children + introns # join lists

                # handle gene components
                for feature in all_gene_features:
                    tracts = [
                        r for r in snp_reader.fetch(
                        feature.chrom, feature.start, feature.end)]
                    if not tracts:
                        continue

                    feature_dict[feature.featuretype] += 1
                    eff_bp_dict[feature.featuretype] += get_eff_bp(
                        bam_reader, chrom, feature.start, feature.end)

            # handle intergenic sequence with a second iteration
            intergenic_tracts = [tract for tract in db.interfeatures(
                [gene for gene in db.features_of_type('gene', limit=region)])
                if tract.end > tract.start]
            for feature in tqdm(intergenic_tracts, desc=f'{chrom} intergenic'):
                tracts = [
                    r for r in snp_reader.fetch(
                    feature.chrom, feature.start, feature.end)]
                if not tracts:
                    continue
                feature_dict[feature.featuretype] += 1
                eff_bp_dict[feature.featuretype] += get_eff_bp(
                    bam_reader, chrom, feature.start, feature.end)

            out_dict = {
                'cross': cross,
                'chrom': chrom, 
                'CDS_count': feature_dict['CDS'],
                'CDS_eff_bp': eff_bp_dict['CDS'],
                'utr5_count': feature_dict['five_prime_UTR'],
                'utr5_eff_bp': eff_bp_dict['five_prime_UTR'],
                'utr3_count': feature_dict['three_prime_UTR'],
                'utr3_eff_bp': eff_bp_dict['three_prime_UTR'],
                'intronic_count': feature_dict['inter_CDS_CDS'],
                'intronic_eff_bp': eff_bp_dict['inter_CDS_CDS'],
                'intergenic_count': feature_dict['inter_gene_gene'],
                'intergenic_eff_bp': eff_bp_dict['inter_gene_gene']}
            writer.writerow(out_dict)
            

def main():
    bam, gff, snp_read_counts, out = args()
    denominator_counts(bam, gff, snp_read_counts, out, lengths)

if __name__ == '__main__':
    main()

        

