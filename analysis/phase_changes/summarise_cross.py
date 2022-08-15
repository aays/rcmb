"""
summarise_cross.py - summarise putative calls into a table
"""

import os
import sys
import csv
import subprocess
import numpy as np
from readcomb.classification import Pair
import argparse
import pysam
from tqdm import tqdm
from multiprocessing import Queue
from multiprocessing import Process

def arg_parser():
    parser = argparse.ArgumentParser(
        description='summarise putative recombinant calls', 
        usage='python summarise_cross.py [options]')

    parser.add_argument('-f', '--bam', required=True,
        type=str, help='BAM output from readcomb-filter and readcomb-fp')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF with calls for cross of interest')
    parser.add_argument('-m', '--mask_size', required=True,
        type=int, help='Masking size')
    parser.add_argument('-p', '--processes', required=False, default=1,
        type=int, help='Number of processes to use (default 1)')
    parser.add_argument('-q', '--mapq', required=False,
        default=50, type=int, help='MAPQ filter')
    parser.add_argument('--base_qual', required=False,
        default=50, type=int, help='base qual filter')
    parser.add_argument('--remove_uninformative', required=False,
        action='store_true', help='Ignore all no phase change reads')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    return parser

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class Processor(Process):
    """
    inherited from multiprocessing.Process

    writes desired pair info to a csv file - each process gets its own temp csv
    file before the files are concatenated

    processes are given pairs of reads to then assemble `rc.Pair` objects from
    and then classify + output quality metrics
    """
    def __init__(self, input_queue, counter_queue, args, p_id, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counter_queue = counter_queue
        self.args = args
        self.p_id = p_id
        self.header = pysam.AlignmentFile(args.bam, 'r').header
        self.fieldnames = [
            'chromosome', 'midpoint', 'start', 'end', 'rel_midpoint', 'call', 'masked_call',
            'mask_size', 'var_count', 'outer_bound', 'min_end_proximity', 'min_vars_in_hap', 
            'var_skew', 'mismatch_var_ratio', 'var_per_hap', 'gc_length', 'indel_count',
            'indel_proximity', 'proximate_indel_length', 'read1_mapq', 'read2_mapq',
            'min_var_qual', 'min_var_depth', 'avg_diff_parental_gq', 'avg_hap_var_gq', 
            'avg_hap_var_depth', 'min_base_qual', 'mean_base_qual', 'read1_length', 
            'read2_length', 'effective_length', 'detection', 'read_name']

    def run(self):

        with open(self.args.out + f'{self.p_id}', 'w', newline='') as f:
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=self.fieldnames)
            writer.writeheader()
            current_pair = self.input_queue.get(block=True)
            while current_pair:
                self.counter_queue.put(1) # increment counter
                record_1 = pysam.AlignedSegment.fromstring(current_pair[0], self.header)
                record_2 = pysam.AlignedSegment.fromstring(current_pair[1], self.header)
                current_pair = self.input_queue.get(block=True)
                if not current_pair:
                    break

                pair = Pair(record_1, record_2, self.args.vcf)

                with HiddenPrints():
                    if all([pair.rec_1.mapq >= self.args.mapq, pair.rec_2.mapq >= self.args.mapq]):
                        pair.classify(masking=self.args.mask_size, quality=self.args.base_qual)
                    else:
                        current_pair = self.input_queue.get(block=True)
                        continue

                if self.args.remove_uninformative and pair.call == 'no_phase_change':
                    current_pair = self.input_queue.get(block=True)
                    continue
                else:
                    start = pair.rec_1.reference_start
                    end = pair.rec_2.reference_start + len(pair.segment_2)
                    if not 'N' in [hap for hap, _, _, _ in pair.detection]:
                        pair_vars = pair.variants_filt
                        detection_indices = [int(hap) - 1 for hap, _, _, _ in pair.detection]
                    else:
                        # remove no match variant from detection and variants_filt
                        no_match_idx = [i for i, var in enumerate(pair.detection) if var[0] == 'N']
                        detection_indices = [
                            int(hap) - 1 for hap in 
                            [hap for hap, _, _, _ in pair.detection if hap != 'N']]
                        pair_vars = [
                            v for i, v in enumerate(pair.variants_filt)
                            if i not in no_match_idx]
                    out_dict = {
                        'chromosome': pair.rec_1.reference_name, 
                        'midpoint': pair.midpoint,
                        'start': pair.rec_1.reference_start,
                        'end': pair.rec_2.reference_end,
                        'rel_midpoint': pair.relative_midpoint,
                        'call': pair.call, 
                        'masked_call': pair.masked_call,
                        'mask_size': self.args.mask_size, 
                        'var_count': len(pair.variants_filt),
                        'outer_bound': pair.outer_bound,
                        'min_end_proximity': pair.min_end_proximity,
                        'min_vars_in_hap': pair.min_variants_in_haplotype,
                        'var_skew': pair.variant_skew, 
                        'mismatch_var_ratio': pair.mismatch_variant_ratio,
                        'var_per_hap': pair.variants_per_haplotype,
                        'gc_length': pair.gene_conversion_len,
                        'indel_count': len(pair.indels) if pair.indels else 0,
                        'indel_proximity': pair.indel_proximity,
                        'proximate_indel_length': pair.proximate_indel_length,
                        'read1_mapq': pair.rec_1.mapq, 'read2_mapq': pair.rec_2.mapq,
                        'min_base_qual': pair.min_base_qual,
                        'mean_base_qual': pair.mean_base_qual,
                        'read1_length': len(pair.rec_1.query_sequence),
                        'read2_length': len(pair.rec_2.query_sequence),
                        'effective_length': end - start,
                        'detection': pair.detection,
                        'read_name': pair.rec_1.query_name
                        }

                    if pair_vars:
                        out_dict['min_var_qual'] = round(min([v.QUAL for v in pair_vars]), 2)
                        out_dict['min_var_depth']: min(np.concatenate([v.gt_depths for v in pair_vars]))
                        out_dict['avg_diff_parental_gq'] = round(
                            sum([abs(v.gt_quals[0] - v.gt_quals[1]) 
                            for v in pair_vars]) / len(pair_vars), 2)
                        out_dict['avg_hap_var_gq'] = round(
                            sum(v.gt_quals[detection_indices[i]]
                            for i, v in enumerate(pair_vars)) / len(pair_vars), 2)
                        out_dict['avg_hap_var_depth'] = sum(v.gt_depths[detection_indices[i]]
                            for i, v in enumerate(pair_vars)) / len(pair_vars)
                    else:
                        out_dict['min_var_qual'] = -1
                        out_dict['min_var_depth'] = -1
                        out_dict['avg_diff_parental_gq'] = -1
                        out_dict['avg_hap_var_gq'] = -1
                        out_dict['avg_hap_var_depth'] = -1
                    writer.writerow(out_dict)

                current_pair = self.input_queue.get(block=True)
        

class Counter(Process):
    """
    counter object for logging if implementation is ever needed

    currently just creates a tqdm progress bar
    """
    def __init__(self, input_queue, args, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.args = args
    
    def run(self):
        ps = subprocess.Popen(['grep', '-v', '@', self.args.bam], stdout=subprocess.PIPE)
        output = subprocess.run(['wc', '-l'], stdin=ps.stdout, stdout=subprocess.PIPE)
        line_count = int(int(output.stdout.decode('utf-8').rstrip('\n')) / 4)
        self.progress = tqdm(total=line_count)

        count = self.input_queue.get(block=True)

        while count:
            self.progress.update(n=1)
            count = self.input_queue.get(block=True)

        self.progress.close()


def parse_reads_parallel(args):
    """parse reads w/ multiprocessing if needed and write to file

    Parameters
    ----------
    args : argparse.ArgumentParser
        command line args

    Returns
    -------
    None
    """

    print('[readcomb] creating processes')
    processes = []
    input_queues = []

    count_input = Queue()
    counter = Counter(count_input, args, daemon=True)
    counter.start()

    for i in range(args.processes):
        input_queue = Queue()
        input_queues.append(input_queue)
        processes.append(
            Processor(input_queue=input_queue, 
                counter_queue=count_input,
                args=args, p_id=i, daemon=True))
        processes[i].start()

    print('[readcomb] processes created')

    bam = pysam.AlignmentFile(args.bam, 'r')
    prev_record = None
    process_idx = 0

    for record in bam:
        if not prev_record:
            prev_record = record
        else:
            if prev_record.query_name != record.query_name:
                raise ValueError(f'read {record.query_name} not paired')

            pair = [prev_record.to_string(), record.to_string()]

            input_queues[process_idx].put(pair)

            if process_idx < args.processes - 1:
                process_idx += 1
            else:
                process_idx = 0
            prev_record = None

    for i, _ in enumerate(processes):
        input_queues[i].put(None)
        input_queues[i].close()
        processes[i].join()

    count_input.put(None)
    count_input.close()
    counter.join()

    print('[readcomb] concatenating temp files')
    with open(args.out, 'w', newline='') as f_out:
        first = True
        for temp_file in [f'{args.out}{i}' for i, _ in enumerate(processes)]:
            with open(temp_file, 'r', newline='') as f_in:
                reader = csv.DictReader(f_in, delimiter='\t')
                if first:
                    writer = csv.DictWriter(f_out, delimiter='\t',
                        fieldnames=reader.fieldnames)
                    writer.writeheader()
                    first = False
                for line in reader:
                    writer.writerow(line)
            os.remove(temp_file)
                    

def main():
    parser = arg_parser()
    parse_reads_parallel(parser.parse_args())

if __name__ == '__main__':
    main()

        

