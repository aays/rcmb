"""
summarise_chrom_eff_bp.py - summarise chrom eff bp, considering SNP deserts

intended to be run as part of callability.smk and not as a standalone script

from workflow:

input:
    'data/callability/tracts/{cross}.read_counts.tsv'
output:
    'data/callability/eff_bp/{cross}.chrom.tsv'

"""

import os
import csv
import pysam
from tqdm import tqdm

def buf_count_newlines_gen(fname):
    # https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-c
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


def summarise_chrom_eff_bp(read_counts_file, tsv_output):
    with open(str(tsv_output), 'w') as f_out:
        fieldnames = [
            'cross', 'chrom', 'eff_bp', 'total_tract_seq']
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        fname = str(read_counts_file)

        with open(fname, 'r') as f_in:
            total_lines = buf_count_newlines_gen(fname)
            reader = csv.DictReader(f_in, delimiter='\t')
            cross = os.path.basename(fname).rstrip('.read_counts.tsv')
            bam_fname = f'data/alignments/bam_prepped/sorted/{cross}.sorted.bam'
            bam_reader = pysam.AlignmentFile(bam_fname, 'rb')
            d = {}

            for line in tqdm(reader, total=total_lines):
                chrom = line['chrom']
                if chrom not in d:
                    d[chrom] = {}
                    d[chrom]['eff_bp'] = 0
                    d[chrom]['total_tract_seq'] = 0

                # ignore SNP deserts
                if int(line['pairs_span']) == 0:
                    continue

                start = int(line['snp_1'])
                end = int(line['snp_2'])
                for p_column in bam_reader.pileup(chrom, start, end, truncate=True):
                    for p_read in p_column.pileups:
                        d[chrom]['eff_bp'] += 1
                d[chrom]['total_tract_seq'] += end - start

            for chrom in sorted(d.keys()):
                writer.writerow({
                    'cross': cross,
                    'chrom': chrom,
                    'eff_bp': d[chrom]['eff_bp'],
                    'total_tract_seq': d[chrom]['total_tract_seq']
                    })

summarise_chrom_eff_bp(snakemake.input, snakemake.output)
