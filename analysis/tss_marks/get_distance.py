"""
get_distance.py - get distance to TSS + chromatin marks
"""

import csv
import math
import pysam
import argparse
from tqdm import tqdm
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='', 
        usage='python script.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='List of recombination events')
    parser.add_argument('-t', '--tss', required=True,
        type=str, help='Path to TSS summary file')
    parser.add_argument('-p', '--peaks', required=True,
        type=str, help='Path to narrowPeak file with chromatin marks')
    # parser.add_argument('--cov_dir', required=True,
        # type=str, help='Dir containing parental coverage files')
    parser.add_argument('-b', '--bam', required=True,
        type=str, help='Path to position-sorted recombinant BAM')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='Path to cross VCF')
    parser.add_argument('-w', '--window_size', required=True,
        type=int, help='Window size')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args

def get_tss_info(chrom, co_pos, tss_tabix, args):
    """
    returns:
        - TSS record
        - dist to TSS
        - intergenic tract size
    """
    tss_chrom = tss_tabix.fetch(chrom)
    tss_positions = []
    tts_positions = []
    for line in tss_chrom:
        chrom, tss_pos, strand, start, end = line.split('\t')
        if strand == '+':
            tss_positions.append(int(start))
            tts_positions.append(int(end))
        elif strand == '-':
            tss_positions.append(int(end))
            tts_positions.append(int(start))
    min_pos_dist = min([abs(pos - co_pos) for pos in tss_positions])
    min_pos = [pos for pos in tss_positions
        if abs(pos - co_pos) == min_pos_dist][0]
    # get intergenic tract size
    intergenic_dists = [(min_pos - tts) for tts in tts_positions]
    # remove any instances where gene's own tts is closer than prev tts
    intergenic_dists = [dist for dist in intergenic_dists if dist > 0]
    # get length
    if len(intergenic_dists) == 0: # first TSS of the chromosome
        intergenic_tract_length = min_pos # 0 -> current TSS
        prev_tts = 0
    else:
        intergenic_tract_length = min([dist for dist in intergenic_dists if dist != 0])
        prev_tts = min_pos - intergenic_tract_length
    
    # get full tss line for later info
    tss_line = [line for line in tss_tabix.fetch(chrom)
        if int(line.split('\t')[1]) == min_pos]

    return min_pos, min_pos_dist, intergenic_tract_length, prev_tts, tss_line

def get_peak_info(chrom, co_pos, peaks_tabix, args):
    """
    returns:
        - start of nearest peak
        - end of nearest peak
    """
    peak_chrom = peaks_tabix.fetch(chrom)
    peak_midpoints_all = []
    for line in peak_chrom:
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        peak_midpoint = int(start + (end - start))
        peak_midpoints_all.append(peak_midpoint)
    min_peak_dist = min([abs(peak_midpoint - co_pos) 
        for peak_midpoint in peak_midpoints_all])
    peak_chrom = peaks_tabix.fetch(chrom)
    for line in peak_chrom:
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        peak_midpoint = int(start + (end - start))
        if abs(peak_midpoint - co_pos) == min_peak_dist:
            nearest_peak = line
            break
    chrom, peak_start, peak_end = nearest_peak.split('\t')[:3]
    return peak_start, peak_end

def get_dist_window(chrom, feature_midpoint, min_pos_dist, args):
    """
    use dist to TSS/peak + feature midpoint + feature pos to get a window and return
        - sum of possible rcmb seq in bp
        - count of SNPs in window

    eg if dist is -125 and windowsize is 100, then
    CO is assigned to the -100 -> -200 bin as a single CO
    out of X possible rcmb seq and based off of Y SNPs 

    calculating X and Y requires getting genomic positions relative to feature
    midpoint

    SNP density is then just Y / window_size,
    and effective denominator is just sum of read pileups, ignoring overlaps
    """
    # get window
    window_start = math.floor(min_pos_dist / args.window_size) * args.window_size
    window_end = window_start + args.window_size

    # get genomic positions of window
    genomic_start = feature_midpoint + window_start
    genomic_end = feature_midpoint + window_end
    if genomic_start < 0:
        print('what the fuck', feature_midpoint, window_start, min_pos_dist, genomic_start)

    # get sum of possible rcmb seq
    bam = pysam.AlignmentFile(args.bam, 'rb')
    eff_bp = 0
    for p_column in bam.pileup(chrom, genomic_start, genomic_end):
        # tested this and it only counts overlapping reads once - nice
        for p_read in p_column.pileups: # all unique read pairs at base
            eff_bp += 1

    # get count of SNPs in window 
    # can be used to get SNP density at whichever scale (per cross, overall, etc)
    # by dividing by sequence considered 
    vcf_reader = VCF(args.vcf)
    snp_count = len([snp for snp in vcf_reader(f'{chrom}:{genomic_start}-{genomic_end}')])
    vcf_reader.close() 

    return window_start, window_end, eff_bp, snp_count


def parse_cos(args):
    with open(args.fname, 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')

        with open(args.out, 'w') as f_out:
            fieldnames = [
                'cross', 'type', 'chromosome', 'midpoint', 'start', 'end',
                'read_name', 'prob', 'tss_window_start', 'tss_window_end', 
                'tss_nearest', 'tss_dist', 'intergenic_tract', 'prev_tts', 
                'tss_snp_count', 'tss_eff_bp_denom', 'peak_dist', 'peak_nearest_start', 
                'peak_nearest_end', 'peak_window_start', 'peak_window_end', 
                'peak_snp_count', 'peak_eff_bp_denom']
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            tss_reader = pysam.TabixFile(args.tss)
            peak_reader = pysam.TabixFile(args.peaks)

            for line in tqdm(reader):
                chrom = line['chromosome']
                if chrom in ['cpDNA', 'mtDNA']:
                    continue
                midpoint = int(line['midpoint'].split(':')[1])

                # tss
                tss_nearest, tss_dist, intergenic_tract_length, prev_tts, tss_line = \
                    get_tss_info(chrom, midpoint, tss_reader, args)
                tss_window_start, tss_window_end, tss_eff_bp, tss_snp_count = \
                    get_dist_window(chrom, tss_nearest, tss_dist, args)

                # peaks
                peak_start, peak_end = get_peak_info(
                    chrom, midpoint, peak_reader, args)
                peak_midpoint = int(peak_start) + int((int(peak_end) - int(peak_start)) / 2)
                peak_dist = int(midpoint - peak_midpoint)
                peak_window_start, peak_window_end, peak_eff_bp, peak_snp_count = \
                    get_dist_window(chrom, peak_midpoint, peak_dist, args)
                out_dict = {
                    'cross': line['cross'], 'type': line['type'],
                    'chromosome': chrom, 'midpoint': line['midpoint'],
                    'start': line['start'], 'end': line['end'],
                    'read_name': line['read_name'], 'prob': line['prob'],
                    'tss_window_start': tss_window_start,
                    'tss_window_end': tss_window_end,
                    'tss_nearest': tss_nearest,
                    'tss_dist': tss_dist,
                    'intergenic_tract': intergenic_tract_length,
                    'prev_tts': prev_tts,
                    'tss_snp_count': tss_snp_count,
                    'tss_eff_bp_denom': tss_eff_bp,
                    'peak_dist': peak_dist,
                    'peak_nearest_start': peak_start,
                    'peak_nearest_end': peak_end,
                    'peak_window_start': peak_window_start,
                    'peak_window_end': peak_window_end,
                    'peak_snp_count': peak_snp_count,
                    'peak_eff_bp_denom': peak_eff_bp
                    }
                writer.writerow(out_dict)

def main():
    parse_cos(args())

if __name__ == '__main__':
    main()

        

