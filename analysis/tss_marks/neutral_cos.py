"""
neutral_cos.py - generate neutral CO expectations 
"""

import csv
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='generate neutral CO expectations', 
        usage='python neutral_cos.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='File containing COs')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.out

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
    'chromosome_17': 6954842,
    'cpDNA': 205535,
    'mtDNA': 15789 # including organelles for mt leakage
    }

def sample_get_neutral(chrom, n_cos, lengths):
    """get evenly spaced CO positions for chrom

    Parameters
    ----------
    chrom : str
        chromosome of interest
    n_cos : int
        number of COs on chromosome
    lengths : dict
        dict containing length of each chromosome


    Returns
    -------
    uniform_sites : list
        list containing positions of evenly spaced sites
    """
    if n_cos == 0:
        return None
    step = int(lengths[chrom] / n_cos)
    uniform_sites = [
        int(step / 2) + (i*step) for i in range(n_cos)]
    return uniform_sites

def generate_neutral(fname, out, lengths):
    """generate evenly spaced COs

    Parameters
    ----------
    fname : str
        path to file containing COs
    out : str
        path to file to write to


    Returns
    -------
    None
    """
    with open(fname, 'r') as f_in:
        with open(out, 'w') as f_out:
            reader = csv.DictReader(f_in, delimiter='\t')
            fieldnames = [
                'cross', 'type', 'chromosome', 'midpoint', 'start', 'end', 
                'read_name', 'prob']
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            
            # get CO counts per chr per sample
            d = {}
            types = {}
            for line in tqdm(reader, desc='co counts'):
                if line['cross'] not in types:
                    types[line['cross']] = line['type']

                if line['chromosome'] in ['cpDNA', 'mtDNA']:
                    continue

                if not line['cross'] in d:
                    # create chrom
                    d[line['cross']] = {
                        k: 0 for k in [f'chromosome_0{i}' for i in range(10)]}
                    for i in range(10, 18):
                        d[line['cross']][f'chromosome_{i}'] = 0

                elif line['cross'] in d:
                    d[line['cross']][line['chromosome']] += 1
            
            for cross in d:
                for chrom in d[cross]:
                    uniform_positions = sample_get_neutral(
                        chrom, d[cross][chrom], lengths)
                    
                    if uniform_positions:
                        for pos in uniform_positions:
                            out_dict = {
                                'cross': cross,
                                'type': types[cross],
                                'chromosome': chrom,
                                'midpoint': f'{chrom}:{pos}',
                                'start': pos - 175, # 350 avg
                                'end': pos + 175,
                                'read_name': 'null',
                                'prob': 1.0
                                }
                            writer.writerow(out_dict)

def main():
    fname, out = args()
    generate_neutral(fname, out, lengths)

if __name__ == '__main__':
    main()

        

