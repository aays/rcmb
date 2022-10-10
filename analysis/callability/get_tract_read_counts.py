"""
get_tract_read_counts.py - get counts of reads per SNP tract

intended to be run as part of callability.smk and not as a standalone script

from workflow:

input:
    tracts = 'data/callability/tracts/{cross}.tracts.tsv',
    bam = 'data/alignments/bam_prepped/sorted/{cross}.sorted.bam
output:
    'data/callability/tracts/{cross}.read_counts.tsv'
"""

import os
import csv
import pysam
from tqdm import tqdm
from readcomb.classification import Pair
from copy import deepcopy
from collections import defaultdict

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def buf_count_newlines_gen(fname):
    # https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python/68385697#68385697
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count

def get_tract_read_counts(tract_file, bam_input, output):
    with open(str(tract_file), 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        total_lines = buf_count_newlines_gen(str(tract_file))

        with open(str(output), 'w') as f_out:
            fieldnames = deepcopy(reader.fieldnames)
            fieldnames.extend(['pairs_span', 'reads_overlap_total'])
            writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()

            bam_reader = pysam.AlignmentFile(str(bam_input))

            for line in tqdm(reader, total=total_lines):
                cross, chrom, snp_1, snp_2, length = line.values()
                snp_1, snp_2 = int(snp_1), int(snp_2)
                out_dict = line

                region = bam_reader.fetch(chrom, snp_1, snp_2)
                reads = [read for read in region]
                if not reads:
                    out_dict['pairs_span'] = 0
                    out_dict['reads_overlap_total'] = 0
                    writer.writerow(out_dict)
                    continue
                else:
                    pairs_span_count = 0
                    reads_overlap_total = 0
                    read_names = [read.query_name for read in reads]

                    # speed up - check if even one read spans
                    min_check = any([
                        read.reference_start <= snp_1 and read.reference_end >= snp_2
                        for read in reads])
                    if min_check:
                        out_dict['pairs_span'] = len(reads)
                        out_dict['reads_overlap_total'] = int(len(reads) / 2) # proxy
                        writer.writerow(out_dict)
                        continue

                    spans = []
                    pairs = read_pair_generator(
                        bam_reader, region_string=f'{chrom}:{snp_1}-{snp_2}')
                    for read1, read2 in pairs:
                        reads_overlap_total += 1
                        spans.append((read1.reference_start, read2.reference_end))
                    for span_start, span_end in spans:
                        if span_start <= snp_1 and span_end >= snp_2:
                            pairs_span_count += 1

                    out_dict['pairs_span'] = pairs_span_count
                    out_dict['reads_overlap_total'] = reads_overlap_total
                    writer.writerow(out_dict)

get_tract_read_counts(
    snakemake.input.tracts, snakemake.input.bam, snakemake.output)


