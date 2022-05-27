"""
summarise_cross.py - summarise putative calls into a table
"""

import os
import sys
import csv
import readcomb.classification as rc
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='summarise putative recombinant calls', 
        usage='python summarise_cross.py [options]')

    parser.add_argument('-f', '--bam', required=True,
        type=str, help='BAM output from readcomb-filter and readcomb-fp')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF with calls for cross of interest')
    parser.add_argument('-m', '--mask_size', required=True,
        type=int, help='Masking size')
    parser.add_argument('--remove_uninformative', required=False,
        action='store_true', help='Ignore all no phase change reads')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.bam, args.vcf, args.mask_size, \
        args.remove_uninformative, args.out

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def parse_reads(bam, vcf, mask_size, remove_uninformative, out):
    """parse reads and write to file

    Parameters
    ----------
    bam : str
        BAM output from readcomb-filter and readcomb-fp
    vcf : str
        path to VCF with calls for cross of interest        
    mask_size : int
        size of region to use for read masking        
    remove_uninformative : bool
        if True, remove all no_phase_change reads
    out : str
        file to write to


    Returns
    -------
    None
    """
    fieldnames = [
        'chromosome', 'midpoint', 'rel_midpoint', 'call', 'masked_call',
        'mask_size', 'var_count', 'min_vars_in_hap', 'var_skew',
        'mismatch_var_ratio', 'var_per_hap', 'gc_length', 'indel_proximity',
        'detection', 'read_name']

    with open(out, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        reader = rc.pairs_creation(bam, vcf)

        for pair in tqdm(reader):
            with HiddenPrints():
                pair.classify(masking=mask_size)

            if remove_uninformative and 'no_phase_change' in [pair.call, pair.masked_call]:
                continue
            else:
                writer.writerow({
                    'chromosome': pair.rec_1.reference_name, 
                    'midpoint': pair.midpoint,
                    'rel_midpoint': pair.relative_midpoint,
                    'call': pair.call, 
                    'masked_call': pair.masked_call,
                    'mask_size': mask_size, 
                    'var_count': len(pair.variants_filt),
                    'min_vars_in_hap': pair.min_variants_in_haplotype,
                    'var_skew': pair.variant_skew, 
                    'mismatch_var_ratio': pair.mismatch_variant_ratio,
                    'var_per_hap': pair.variants_per_haplotype,
                    'gc_length': pair.gene_conversion_len,
                    'indel_proximity': pair.indel_proximity,
                    'detection': pair.detection,
                    'read_name': pair.rec_1.query_name
                    })

def main():
    parse_reads(*args())

if __name__ == '__main__':
    main()

        

