"""
biased_seg.py - get counts of reads matching either ancestor across windows
"""

import csv
import argparse
import math
import numpy as np
import pysam
import readcomb.classification as rc
from multiprocessing import Queue
from multiprocessing import Process
from tqdm import tqdm

def arg_parser():
    parser = argparse.ArgumentParser(
        description='get counts of reads matching either ancestor', 
        usage='python biased_seg.py [options]')

    parser.add_argument('-f', '--bam', required=True,
        type=str, help='Read name sorted sam/bam file containing reads of interest')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF with calls for cross of interest')
    parser.add_argument('-w', '--window_size', required=True, default=1000,
        type=int, help='Windowsize (default 1000)')
    parser.add_argument('-q', '--base_qual', required=False, default=0,
        type=int, help='Base qual filter (default 0)')
    parser.add_argument('-p', '--processes', required=False, default=1,
        type=int, help='Number of processes (default 1)')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    return parser

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

def create_array_set(lengths, window_size):
    """generate generic set of arrays for provided chromosomes/contigs

    given a window size, will create length-L numpy arrays for each provided chromosome
    where L is chrom_length / window_size. the individual elements of the array
    represent windows, and windows can be recapitulated via element_index * window_size

    Parameters
    ----------
    lengths : dict
        dict of chromosome/contig/scaffold lengths
    window_size : int
        desired window size to split chrs/contigs/scaffolds into

    Returns
    -------
    arrays : dict
        dict of numpy arrays, with keys corresponding to chrs/contigs/scaffolds

    """
    arrays = {}
    for chrom in lengths.keys():
        array_len = len(range(0, lengths[chrom], window_size)) + 1 # add 1 for window end
        arrays[chrom] = np.zeros(array_len)
    return arrays


class Processor(Process):
    def __init__(self, input_queue, counter_queue, lengths, args, p_id, **kwargs):
        """
        inherited from multiprocessing.Process

        spawns 2n arrays, where n = number of chromosomes (hardcoded at 17 right
        now, might add an option to modify later) - one set per parent

        chromosome lengths are also hardcoded for now but could easily be
        updated as an input arg in the future

        iterates through read pairs passed through it and assigns ancestry - will
        then increment corresponding window for corresponding parent
        """
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.counter_queue = counter_queue
        self.output_queue = output_queue
        self.args = args
        self.p_id = p_id # TODO: not sure I need this here
        self.header = pysam.AlignmentFile(args.bam, 'r').header
        self.lengths = lengths
        self.parent1_arrays = create_array_set(lengths)
        self.parent2_arrays = create_array_set(lengths)
        self.removed_arrays = create_array_set(lengths)
        self.failed_mapq = create_array_set(lengths)

    def run(self):
        current_pair = self.input_queue.get(block=True)
        while current_pair:
            self.counter_queue.put(1) # increment counter
            record_1 = pysam.AlignedSegment.fromstring(current_pair[0], self.header)
            record_2 = pysam.AlignedSegment.fromstring(current_pair[1], self.header)

            pair = Pair(record_1, record_2, self.args.vcf)
            pair.classify(masking=0, quality=self.args.base_qual)
            chrom = pair.rec_1.reference_name
            window = int(pair.midpoint.split(':')[1])

            # mapq filter if provided
            if all([pair.rec_1.mapq < self.args.mapq, pair.rec_2.mapq < self.args.mapq]):
                self.failed_mapq[chrom][window] += 1
                current_pair = self.input_queue.get(block=True)
                continue

            # ignore recombinant reads, reads with no variants, and reads on contigs
            if (
                pair.call != 'no_phase_change' or 
                len(pair.detection) == 0 or
                pair.rec_1.query_name not in lengths.keys()
            ):
                self.removed_arrays[chrom][window] += 1
                current_pair = self.input_queue.get(block=True)
                continue

            elif pair.call == 'no_phase_change':
                hap = list(set([hap for hap, _, _, _ in pair.detection if hap != 'N']))[0]
                assert len(hap) == 1
                window = math.floor(int(pair.midpoint.split(':')[1]) / self.args.window_size)
                if hap == '1':
                    self.parent1_arrays[record_1.reference_name][window] += 1
                elif hap == '2':
                    self.parent2_arrays[record_1.reference_name][window] += 1
                current_pair = self.input_queue.get(block=True)
                continue

            else:
                raise Exception("this probably shouldn't happen")

        self.output_queue.put(
            [self.parent1_arrays, self.parent2_arrays, self.removed_arrays, self.failed_mapq])

class Counter(Process):
    """
    counter object for logging if implementation is ever needed
    currently just creates a tqdm progress bar
    """
    def __init__(self, input_queue, args, lengths, **kwargs):
        Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.args = args
        self.lengths = lengths
    
    def run(self):
        stats = subprocess.check_output(['samtools', 'idxstats', args.bam])
        stats_table = pd.read_csv(
            BytesIO(stats), sep='\t',
            names=['chrom', 'length', 'map_reads', 'unmap_reads'])
        stats_table = stats_table[stats_table.chrom.isin(lengths.keys())]
        pairs_sum = stats_table['map_reads'].sum()
        pair_count = int(pairs_sum / 2) # divide 2 to get pairs

        self.progress = tqdm(total=pair_count)
        count = self.input_queue.get(block=True)

        while count:
            self.progress.update(n=1)
            count = self.input_queue.get(block=True)

        self.progress.close()


def biased_seg_parallel(args, lengths):
    """calculate biased seg stats using multiprocessing

    function will spawn n processes and distribute reads to each -
    after all reads have been iterated over, will concatenate resulting
    arrays and write to file

    Parameters
    ----------
    args : argparse.ArgumentParser
        command line args
    lengths : dict
        dictionary containing chromosome lengths


    Returns
    -------
    None
    """

    print('[rcmb] creating processes')
    processes = []
    input_queues = []
    output_queues = []

    count_input = Queue()
    counter = Counter(count_input, args, lengths, daemon=True)

    for i in range(args.processes):
        input_queue = Queue()
        input_queues.append(input_queue)
        output_queue = Queue()
        processes.append(
            Processor(input_queue=input_queue,
                output_queue=output_queue,
                counter_queue=count_input,
                args=args, p_id=i, daemon=True))
        processes[i].start()

    print(f'[rcmb] {len(processes)} processes created')

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

    # close processes
    for i, _ in enumerate(processes):
        input_queues[i].put(None)
        input_queues[i].close()
        processes[i].join()

    count_input.put(None)
    count_input.close()
    counter.join()

    # merge arrays
    def concat_arrays(main_array_set, current_array_set):
        out_array_set = {}
        # elementwise addition for each chrom array
        for chrom in main_array_set:
            out_array_set[chrom] = main_array_set[chrom] + current_array_set[chrom]
        return out_array_set

    parent1_all = create_array_set(lengths)
    parent2_all = create_array_set(lengths)
    removed_all = create_array_set(lengths)
    failed_mapq_all = create_array_set(lengths)

    for out_queue in output_queues:
        parent1_arrays, parent2_arrays, removed, failed_mapq = out_queue.get()
        parent1_all = concat_arrays(parent1_all, parent1_arrays)
        parent2_all = concat_arrays(parent2_all, parent2_arrays)
        removed_all = concat_arrays(removed_all, removed)
        failed_mapq_all = concat_arrays(failed_mapq_all, failed_mapq)
        
    # write to file
    with open(args.out, 'w') as f:
        fieldnames = [
            'chromosome', 'window_start', 'window_end', 'parent1_pairs', 'parent2_pairs', 
            'total_kept_pairs', 'parent1_perc', 'parent2_perc', 'pairs_skipped', 'failed_mapq']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for chrom in lengths:
            windows = range(0, lengths[chrom], args.window_size)
            for i, window_start in enumerate(windows):
                window_end = window_start + args.window_size \
                    if i < len(windows) - 1 else lengths[chrom]

                out_dict = {
                    'chromosome': chrom,
                    'window_start': window_start,
                    'window_end': window_end,
                    'parent1_pairs': parent1_all[chrom][i],
                    'parent2_pairs': parent2_all[chrom][i],
                    'pairs_skipped': removed_all[chrom][i],
                    'failed_mapq': failed_mapq_all[chrom][i]
                    }

                out_dict['total_kept_pairs'] = parent1_all[chrom][i] + parent2_all[chrom][i]
                out_dict['parent1_perc'] = out_dict['parent1_pairs'] / out_dict['total_kept_pairs']
                out_dict['parent2_perc'] = out_dict['parent2_pairs'] / out_dict['total_kept_pairs']

                writer.writerow(out_dict)

def main():
    parser = arg_parser()
    args = parser.parse_args()
    biased_seg_parallel(args, lengths)

if __name__ == '__main__':
    main()

        

