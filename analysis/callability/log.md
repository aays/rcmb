
## 19/9/2022

today - need to estimate callability for each cross

two components:

1. get inter-SNP tract lengths - those above ~350 are likely not callable
2. get coverage distributions for each parent, since anomalies could mean deletions

let's start with inter-SNP tracts - for each VCF, just need to get the
distance from each SNP to the next one and report it

I could do this by generating some sort of bed file along the lines of

```
cross chrom snp_1 snp_2 dist
```

for each individual cross - let's get a script going that does just that,
using `itertools.pairwise` - although that's not implemented yet
in the version of python I have, so:

```python
# from itertools documentation
import itertools
def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

# seems to work as intended...
from cyvcf2 import VCF
reader = VCF('data/genotyping/vcf_filtered/2343x1952.vcf.gz')
p = pairwise(reader)
l = []
l.append(next(p))
l.append(next(p))
l
# [(Variant(chromosome_01:18618 T/G), Variant(chromosome_01:19757 A/C)), 
# (Variant(chromosome_01:19757 A/C), Variant(chromosome_01:23577 C/T))]
```

going to create a quick script called `create_snp_tract_bed.py`
for use with a snakemake workflow -

```bash
time snakemake -pr -s analysis/callability/callability.smk --cores 12
# completed in ~20 min
```

## 21/9/2022

today - just calculating the 'fudge factors' from these tracts

going to summarise tracts with a file containing counts, like so:

```
cross tract_length count
2343x1691 100 24
2343x1691 101 73
2343x1691 102 105
```

and so on

let's update the workflow real quick:

```bash
time snakemake -pr -s analysis/callability/callability.smk --cores 12
```

next up - looking for coverage abnormalities in the parents

going to take a pretty simple approach to this - will slide a 2 kbp window
by 1 kbp each time and get the coverage of it

for each window, need to calc -

1. number of reads
2. number of reads in left window
3. number of reads in right window
4. ratio of reads vs average # of reads in flanking windows

I could use samtools coverage to do the 'hard part' and then slide
across windows to calculate rate ratios etc:

```python
import csv
import subprocess

lengths = {...} # get chrom lengths

fname = ...

fieldnames = [
    '#rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage',
    'meandepth', 'meanbaseq', 'meanmapq']

cmd = 'samtools coverage -H -r {chrom}:{start}-{end} {fname}'

with open('out.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(fieldnames)
    for chrom in lengths:
        for window_start in range(0, lengths[chrom], 1000):
            window_end = window_start + 2000
            if window_end > lengths[chrom]:
                window_end = lengths[chrom]
            proc = subprocess.run(
                cmd.format(
                    chrom=chrom,
                    start=window_start,
                    end=window_end,
                    fname=fname).split(' '),
                check=True,
                stdout=subprocess.PIPE)
        writer.writerow(proc.stdout.decode('utf-8').rstrip('\n').split('\t'))

```

going to generalize this to a snakemake script and then get to it:

```bash
# will run windowed_coverage rule
time snakemake -pr -s analysis/callability/callability.smk --cores 16
```

## 8/10/2022

I just realized there's a far better way to do this - for each of the tracts,
I should retrieve the number of reads that cover it (eg span from before left bound
to past right bound)

I can just query the tract files in `data/callability/tracts` directly and create
`.counts.tsv` files - just need a simple script that will, for each cross,
fetch all reads at the SNP bound coordinates and count the ones that span

giving this a go:

```bash
time snakemake -pr -s analysis/callability/callability.smk --cores 16
# took 4 hours
```

## 9/10/2022

today - using these files to get correct rates 

going to generate callable denominators by looking at tracts > 0 -
let's just generate files containing per-chrom counts, since these can be summed afterwards too

going to add this as a snakemake script 

```bash
time snakemake -pr -s analysis/callability.smk
# run summarise_tracts rule
```

I need to also get per-chr effective bp - this needs to be
taking into account SNP deserts, so I need to scan the genome
and not count effective bp that happens to be in those

basically - for each tract that has more than one read spanning,
use pileup with the start/end coord to get the total amount of bp
that's actually within those tracts specifically

see `tss_marks/get_distance.py` for the pileup thing 

here goes:

```bash
time snakemake -pr -s analysis/callability.smk --cores 16
# runs summarise_chrom_eff_bp rule
```





