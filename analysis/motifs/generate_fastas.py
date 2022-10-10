"""
generate_fastas.py - generate FASTA files containing recombinant sequence
"""

import csv
import argparse
import numpy as np
from tqdm import tqdm
from Bio import SeqIO

def args():
    parser = argparse.ArgumentParser(
        description='generate fasta files w/ recombinant sequence', 
        usage='python generate_fastas.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='Crossovers file')
    parser.add_argument('-w', '--window_size', required=True,
        type=int, help='Amount of sequence to retrieve on either side of cross midpoint')
    parser.add_argument('-d', '--fasta_dir', required=True,
        type=str, help='Dir containing sample-chrom-separated FASTAs')
    parser.add_argument('-r', '--random', required=False,
        action='store_true', help='Generate equivalent amount of random sequences')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to (.fa/.fasta)')

    args = parser.parse_args()

    return args.fname, args.window_size, args.fasta_dir, args.random, args.out

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

def get_seq(parent1, parent2, line, window_size, fasta_dir, random, lengths):
    """docstring

    Parameters
    ----------
    parent1, parent2 : str
        strings containing parental strains
    line : dict
        line output from csv.DictReader
    window_size : int
        window size in bp of sequence to include on either side
    fasta_dir : str
        path to dir containing sample-chr fastas
    random : bool
        if True, return randomly drawn seq from same chr (for null expectation)
    lengths : dict
        dict of chrom lengths

    Returns
    -------
    total_seq : str
        combined seq of length window_size * 2
    midpoint : int
        midpoint coord of sequence
    """
    if not parent1.startswith('GB') and not parent1.startswith('CC'):
        parent1 = 'CC' + parent1
    if not parent2.startswith('GB') and not parent2.startswith('CC'):
        parent2 = 'CC' + parent2
    chrom = line['chromosome']
    if not random:
        midpoint = int(line['midpoint'].split(':')[1])
    elif random:
        # draw random value across chrom
        midpoint_range = [window_size, lengths[chrom] - window_size] # can't be edge of chrom
        midpoint = np.random.randint(midpoint_range[0], midpoint_range[1])
        
    haps = list(set([hap for hap, pos, base, qual in eval(line['detection']) if hap != 'N']))
    if len(haps) > 2:
        print(line)
        raise Exception('wtf')

    # get order of fastas from which to get sequence
    fastas = []
    for hap in haps:
        if hap == '1':
            fastas.append(fasta_dir + f'{parent1}.{chrom}.fa')
        elif hap == '2':
            fastas.append(fasta_dir + f'{parent2}.{chrom}.fa')

    # get positions on either end of midpoint
    coords = [int(midpoint - window_size), midpoint, int(midpoint + window_size)]
    seq_left = str(SeqIO.read(fastas[0], 'fasta')[coords[0]:coords[1]].seq)
    seq_right = str(SeqIO.read(fastas[1], 'fasta')[coords[1]:coords[2]].seq)
    total_seq = seq_left + seq_right
    return total_seq, midpoint
            

def generate_fastas(fname, window_size, fasta_dir, random, lengths, out):
    """docstring

    Parameters
    ----------
    fname : str
        path to CO list
    window_size : int
        window size in bp of sequence to include on either side
    fasta_dir : str
        path to dir containing sample-chr fastas
    out : str
        file to write to


    Returns
    -------
    None
    """
    if not fasta_dir.endswith('/'):
        fasta_dir += '/'
    with open(out, 'w') as f_out:
        with open(fname, 'r') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')

            for line in tqdm(reader):
                cross = line['cross']
                chrom = line['chromosome']
                read_name = line['read_name']
                if chrom in ['cpDNA', 'mtDNA']:
                    continue
                if cross in ['3071x2931', '3071x3062']:
                    continue
                parent1, parent2 = cross.split('x')
                seq, midpoint = get_seq(
                    parent1, parent2, line, window_size, fasta_dir, random, lengths)
                if not random:
                    f_out.write(f'>{cross}-{chrom}:{midpoint}-{read_name}\n')
                elif random:
                    f_out.write(f'>{cross}-{chrom}:{midpoint}-random\n')
                f_out.write(seq + '\n')

def main():
    fname, window_size, fasta_dir, random, out = args()
    generate_fastas(fname, window_size, fasta_dir, random, lengths, out)

if __name__ == '__main__':
    main()

        

