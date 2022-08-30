
## 5/6/2022

today - plan out segregation distortion landscape analysis

I think this should be done in a windowed fashion across, say, 2 kb windows - 
this is somewhat arbitrary but I figure I could use it as a starting point

I could use `classification.Pair.get_midpoint()` for the midpoints of reads,
and bin them into 2 kb windows that way - this would be a concern since
midpoint calc is dependent on phase change status, but we should actually
stick to 

1. no phase change reads
2. reads with any variants at all (obviously)

since we can't really use phase change reads to look at biased segregation at all, for
obvious reasons

now the main issue here is that reads that can be parsed by readcomb are _not_ actually
sorted by position, which would mean creating some sort of lookup structure or, say, populating
some kind of nested dictionary structure that encompasses windows

I _could_ also 'break' readcomb a bit and manually make 'pairs' out of a sorted BAM
simply by setting rec 1 and rec 2 to the exact same read - this might actually be bonkers enough
to work - although I lose resolution by treating these as single end reads and I don't think
I'd like that very much

I think the nested dictionary may indeed be the play - and I wonder if I could do something
gnarly with temp files that could then be sorted afterwards in some way - need to think
on this more tomorrow

## 12/8/2022

so this was on the back burner for a while, but going to finally get on it now

I've had a while to think about this and I think parallelizing is actually the way,
instead of some lookup stuff - if I parallelize with `multiprocessing`, I can just use
the base `readcomb.classification.Pair` methods/attributes to do the hard part

going to start on a script called `biased_seg.py` (creative!) that will, for each
of n processes (and given a windowsize):

- spawn 34 arrays, one per chrom per parent - size of array is dependent on windowsize
    - eg if windowsize is 1000, array is `chrom_length` / 1000
- classify a `Pair` (TODO: how does `Pair` handle reads with no variants at all? double check)
- if pair has informative variants `pair.call` is `no_phase_change`, get ancestry from `detection`
- increment corresponding window value in array by 1
- once process is done, write arrays to files like so:

```
chrom   window_start    window_end  parent1_pairs   parent2_pairs   parent1_frac
chromosome_01   4000    5000    675 684 0.496
chromosome_01   5000    6000    500 400 0.555
```

and then concatenate somehow 

OR - something like

```
# after processes are done
main_array_collection = {} # 34 empty arrays

for process in processes:
    # get final output array collection
    # increment main array collection

# write main array to file
```

this may be fairly memory hungry - but this way, after the processes are done,
only 34 (from process) + 34 (main) are being held in memory at a time,
and it seems the memory footprint isn't huge:

```python
import numpy as np

x = np.zeros(8000) # chr1 is 8m bp - windows of 1kbp mean 8000 windows
sys.getsizeof(x) # 64096
x.nbytes # 64000

y = np.full(8000, 327489032747) # random large number
sys.getsizeof(y) # 64096
y.nbytes # 64000

# modifying values doesn't impact size
y[0] = 42
sys.getsizeof(y) # 64096
y.nbytes # 64000
```

I think a concatenation operation would probably look like this anyways - so
might as well work with this

## 26/8/2022

alright - spent all of yesterday finalizing this script, and now to give it a go

if I get this working the next step will probably to snakemake-ify it - but one step
at a time

```bash
time python analysis/biased_seg/biased_seg.py \
--bam data/alignments/bam_prepped/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--window_size 1000 \
--base_qual 20 \
--mapq 1 \
--processes 12 \
--out biased_seg_test.tsv
```

after some debugging this looks good! but it's S L O W - tqdm was predicting ~55 hours with just 1
process

probably going to crank it up to ~30 for the full runs - but first let's just do a quick run 
on the smallest file I have

```bash
nice -n 19 time python analysis/biased_seg/biased_seg.py \
--bam data/alignments/bam_prepped/3086x3059.sorted.bam \
--vcf data/genotyping/vcf_filtered/3086x3059.vcf.gz \
--window_size 1000 \
--base_qual 20 \
--mapq 1 \
--processes 20 \
--out biased_seg_test.tsv
```

so, after a few hours of this potentially breaking the server and taking up way too much memory,
it turns out I need to find a way to limit queues - or perhaps spawn new input queues/processes
after a queue has been cleared? 

thought for tomorrow - try looping through the bam, spawning a new process to handle n pairs,
spawning a new process to handle the next n, etc

once a queue for a given process is empty (maybe with the `task_done` thing?) then 
merge the array values from the output queue, kill the process, and create a new one for the
next n? 

I feel like this should be a 'solved' problem - but I guess I need to keep reading

## 29/8/2022

switching it up - going to be doing this with `multiprocessing.Pool` instead since
that combines the queue and process stages of this


