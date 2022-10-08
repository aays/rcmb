
## 26/10/2021

today - getting `readcomb-filter` going with a snakemake workflow, but providing
12 threads

(see alignment log for test run)

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk \
--cores 12 # run rule readcomb_filter
```

## 2/11/2021

this ran fine and completed in 35 hours - log also looks good - but this needs
to be redone after curating the SNP dataset some more

back to the `genotyping` log for now

## 8/11/2021

today - converting these first pass sams to bams and then pulling out some
high quality reads to look at in IGV

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk \
--cores 12 # run rules bam_convert and bam_idx
```

using readcomb on this and there are a few bugs to fix/things to add:

- the `no_match` check for a given `Pair` object tests to see if there are Ns in `detection_1/2`,
but these are lists of tuples - we should be checking for Ns in the tuples, not in the list as a whole
- it would be useful to have some list of tuples linking alleles to parental assignments, for reference
when looking back at the cyvcf records
    - could update `downstream_phase_detection` to incl. allele as third element in tuples
- if both haplotypes only have one supporting variant, `min_variants_in_haplotype` remains `None`
    - is this always the case? need to run through a file and see how often this happens
- it would be useful to have convenience methods that return both `detection_1` and `detection_2` combined -
same for `variants_1` and `variants_2` - could have an optional arg that adds a None in the list
at the boundary of the two reads
- the `variants` attributes still contain variants that are ultimately discarded for being uninformative
- include location attribute that returns full segment spanned by read pair in samtools format
- update midpoint attribute to also return in samtools format

also need to regen the VCF index files in `genotyping/vcf` - keep raising warnings about
the idx files being older, which makes no sense because they should have been generated immediately
after, but it can't hurt to update them regardless

ALSO need to rerun readcomb on a few parental bams - could do something like running on GB119 bam
with GB119x3062 VCF, since the ideal expectation would be

1. no phase changes
2. all variants are GB119 variants (e.g. marked as parent 1)

## 15/11/2021

alright - been a few days of testing now but most of the above things are implemented
and I've fixed a few bugs

next up - trying to rerun on a parental bam - here's a run on a single bam:

```bash
mkdir -p data/phase_changes/parental/

time readcomb-bamprep \
--bam data/alignments/parental_bam/GB119.bam \
--threads 20 \
--outdir data/alignments/bam_filtered/

time readcomb-filter \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119x3062.vcf.gz \
--processes 16 \
--mode phase_change \
--log data/phase_changes/parental.log \
--quality 30 \
--out data/phase_changes/parental/GB119.sam
```

what are these??

```
833040 phase change-containing read pairs from total 9999546 read pairs
2510732 reads had no-match variants
5787913 reads did not have enough variants (> 0) to call
time taken: 0:18:31
```

for one, looks like we have some strange unphased calls that need to be accounted for:

```python
>>> pairs[0].detection
[('2', 1790716, 'C|C'), ('2', 1790718, 'C|C'), ('2', 1790728, 'T|T'), ('2', 1790732, 'T|T'), ('2', 1790733, 'T|T'), ('2', 1790734, 'C|C'), ('2', 1790736, 'A|A'), ('2', 1790752, 'A|A'), ('2', 1790756, 'A|A'), ('2', 1790759, 'T|T'), ('1', 1790777, 'G'), ('N', 1790784, None), ('N', 1790789, None), ('2', 1790936, 'A|A'), ('2', 1790960, 'C|C'), ('2', 1790962, 'G'), ('2', 1790965, 'G|G'), ('2', 1790967, 'G'), ('N', 1790973, None), ('2', 1790977, 'A|A'), ('N', 1790979, None), ('2', 1790983, 'G|G'), ('N', 1790988, None), ('2', 1790992, 'G|G'), ('2', 1791004, 'T'), ('2', 1791023, 'T|T')]

>>> x = pairs[0]
>>> x.variants_filt
[Variant(chromosome_16:1790717 T/C), Variant(chromosome_16:1790719 A/C), Variant(chromosome_16:1790729 C/T), Variant(chromosome_16:1790733 G/T), Variant(chromosome_16:1790734 G/T), Variant(chromosome_16:1790735 A/C), Variant(chromosome_16:1790737 G/A), Variant(chromosome_16:1790753 G/A), Variant(chromosome_16:1790757 G/A), Variant(chromosome_16:1790760 G/T), Variant(chromosome_16:1790778 A/G), Variant(chromosome_16:1790785 T/C), Variant(chromosome_16:1790790 A/G), Variant(chromosome_16:1790937 G/A), Variant(chromosome_16:1790961 A/C), Variant(chromosome_16:1790963 G/GGGGATGCGA), Variant(chromosome_16:1790966 A/G), Variant(chromosome_16:1790968 G/GGA), Variant(chromosome_16:1790974 GAGC/G), Variant(chromosome_16:1790978 G/A), Variant(chromosome_16:1790980 G/A), Variant(chromosome_16:1790984 A/G), Variant(chromosome_16:1790989 C/T,G), Variant(chromosome_16:1790993 A/G), Variant(chromosome_16:1791005 T/C), Variant(chromosome_16:1791024 G/T)]
>>> x.variants_filt[0].gt_bases
array(['C|C', 'T/T'], dtype='<U3')
```

## 2/1/2022

alright, finally getting back to this now that the salt results are done

looks like I have some uncommitted changes where I tried to have readcomb correctly
handle instances of `|` phase separators (as opposed to `/`) - this happens in `classification`
but probably needs to be added to `filter` as well

that said, vcfprep _should_ be removing hets anyways - I've updated the 'splitting' of the
genotypes to account for both `/` and `|` now, so time to give this another go

```bash
time ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119x3062.vcf.gz \
--processes 16 \
--mode phase_change \
--log data/phase_changes/parental.log \
--quality 30 \
--out data/phase_changes/parental/GB119.sam
```

alright - I let this run for about 5 min and it returned thousands of reads already - 
let's try and troubleshoot, given we should have basically none

```bash
head -n 300 data/phase_changes/parental/GB119.sam > data/phase_changes/parental/GB119_filt.sam
```

and then from the readcomb folder so that I can use the updated classification script:

```python
import classification as rc
import pysam

bam_fname = '../../rcmb/data/phase_changes/parental/GB119_filt.sam'
vcf_fname = '../../rcmb/data/genotyping/vcf_filtered/GB119x3062.vcf.gz'

pairs = [pair for pair in rc.pairs_creation(bam_fname, vcf_fname)]
```

seems this third record is matching the other sample:

```python
>>> x.classify()
>>> print(x)
Record name: K00133:299:HHLF2BBXX:5:1101:991:24138#
Read1: chromosome_15:1012883-1012964
Read2: chromosome_15:1013101-1013149
VCF: ../../rcmb/data/genotyping/vcf_filtered/GB119x3062.vcf.gz
Unmatched Variant(s): False
Condensed: [['GB119', 1012883, 1013011], ['CC3062', 1013011, 1013149]]
Call: cross_over
Condensed Masked: [['GB119', 1012953, 1013011], ['CC3062', 1013011, 1013079]]
Masked Call: cross_over
Midpoint: chromosome_15:1013011
Variants Per Haplotype: 2.0
Gene Conversion Length: N/A
>>> x.variants_filt
[Variant(chromosome_15:1012939 C/CTGGTCGCG), Variant(chromosome_15:1013130 A/G), Variant(chromosome_15:1013140 GATGTAACA/G)]
>>> x.variants_filt[2]
Variant(chromosome_15:1013140 GATGTAACA/G)
```

trying to manually retrace the steps in `downstream_phase_detection` for this variant

I figured it out!

```python
>>> x.detection
[('2', 1012938, 'C'), ('2', 1013129, 'G'), ('1', 1013139, 'GATGTAACA')]

>>> x.rec_2.reference_start
1013101

>>> cigar_tuples = x.rec_2.cigartuples
>>> cigar_tuples
[(0, 48)]

>>> x.rec_2.query_qualities.tolist()
[41, 37, 37, 41, 37, 41, 37, 41, 37, 27, 41, 41, 41, 41, 41, 41, 37, 37, 32, 12, 41, 41, 41, 41, 41, 41, 41, 37, 41, 32, 41, 41, 41, 37, 22, 41, 41, 41, 41, 37, 37, 37, 37, 27, 32, 32, 32, 32]

>>> x.rec_2.query_qualities.tolist().index(12)
19

>>> x.variants_filt[2].gt_quals
array([96., 30.], dtype=float32)

>>> x.rec_2.reference_start
1013101

>>> x.rec_2.query_sequence[19:]
'ACTGTGGTCGCGTAACCGTGATGTAACAT' # it's not the base with qual < 30

>>> x.rec_2.query_sequence
'GTGCCATGGATCTGACAACACTGTGGTCGCGTAACCGTGATGTAACAT'

>>> x.rec_2.query_sequence[38:]
'GATGTAACAT'

>>> x.variants_filt
[Variant(chromosome_15:1012939 C/CTGGTCGCG), Variant(chromosome_15:1013130 A/G), Variant(chromosome_15:1013140 GATGTAACA/G)]
```

basically - this is a 0, 48 tuple, which means it's an exact match with the
reference - for some reason, the 3062 call there is an insertion that _still
matches the reference_ (perhaps a small tandem duplication?) but technically
both alleles are 'correct' based on what we see in the read but readcomb is
currently (and I think incorrectly) set to favour the longer allele in the case
of an indel 

I think - when we have this case where both alleles of an indel technically match,
we treat the site as uninformative entirely! 

let's see if this is the case with the next record as well:

```python
>>> x = pairs[1]
>>> x.classify()
>>> x
<classification.Pair object at 0x7f86fc165c18>
>>> print(x)
Record name: K00133:299:HHLF2BBXX:5:1101:991:24947#
Read1: chromosome_04:4051890-4051990
Read2: chromosome_04:4052081-4052181
VCF: ../../rcmb/data/genotyping/vcf_filtered/GB119x3062.vcf.gz
Unmatched Variant(s): False
Condensed: [['GB119', 4051890, 4051912], ['CC3062', 4051912, 4051926], ['GB119', 4051926, 4052181]]
Call: gene_conversion
Condensed Masked: [['GB119', 4051960, 4052111]]
Masked Call: no_phase_change
Midpoint: chromosome_04:4051919
Variants Per Haplotype: 14.0
Gene Conversion Length: 14
>>> x.variants_filt
[Variant(chromosome_04:4051934 G/A), Variant(chromosome_04:4051936 A/AT), Variant(chromosome_04:4051942 G/A), Variant(chromosome_04:4051945 T/C), Variant(chromosome_04:4051958 G/C), Variant(chromosome_04:4051962 T/G), Variant(chromosome_04:4051963 T/G), Variant(chromosome_04:4051964 A/G), Variant(chromosome_04:4052083 C/T), Variant(chromosome_04:4052098 C/A), Variant(chromosome_04:4052140 G/GTGCCCAGCA), Variant(chromosome_04:4052175 C/T), Variant(chromosome_04:4052180 A/T), Variant(chromosome_04:4052181 T/G)]
>>> x.detection
[('2', 4051933, 'G'), ('1', 4051935, 'AT'), ('2', 4051941, 'G'), ('2', 4051944, 'T'), ('2', 4051957, 'G'), ('2', 4051961, 'T'), ('2', 4051962, 'T'), ('2', 4051963, 'A'), ('2', 4052082, 'C'), ('2', 4052097, 'C'), ('2', 4052139, 'G'), ('2', 4052174, 'C'), ('2', 4052179, 'A'), ('2', 4052180, 'T')]

>>> x.variants_filt[1]
Variant(chromosome_04:4051936 A/AT)
>>> x.rec_1.cigartuples
[(0, 100)]
```

same issue! even though this is an exact match with the reference, the longer allele
is favoured (AT over A)

also super odd to see how long the variants per haplotype score is - looks like it's because
this is calculated _with_ masking - which is fine but good to keep in mind (and maybe worth making
more explicit)

could help to also print min variant in haplotype as part of the str representation - implementing
all this now

alright - so takeaways:

- update `vcfprep` so that if there's an indel where the shorter allele is a subset of the longer one,
it's no longer considered an informative allele (going to test the impact of this first by updating
`phase_detection` in `filter` and `classification`
- update `variants_per_haplotype` to make clear that it refers to post-masking
- add `min_variants_in_haplotype` to the string representation

before doing these, just quickly updating `phase_detection` to ignore any indels
where the alleles 'overlap' just to see how this impacts filter's output

```bash
time ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119x3062.vcf.gz \
--processes 16 \
--mode phase_change \
--log data/phase_changes/parental.log \
--quality 30 \
--out data/phase_changes/parental/GB119.sam
```

still thousands turning up... here we go again:

```python
# after head -n 300'ing the above file

import classification as rc
import pysam

bam_fname = '../../rcmb/data/phase_changes/parental/GB119_filt.sam'
vcf_fname = '../../rcmb/data/genotyping/vcf_filtered/GB119x3062.vcf.gz'

pairs = [pair for pair in rc.pairs_creation(bam_fname, vcf_fname)]
```

the culprit still seems to be those indels - I'm just going to remove them
altogether using vcfprep

trying this out:

```bash
zcat data/genotyping/vcf_filtered/GB119x3062.vcf.gz | wc -l # 1656201

time ../readcomb/readcomb/vcfprep.py --vcf data/genotyping/vcf/GB119x3062.vcf.gz \
--no_hets --min_GQ 30 --out data/genotyping/vcf_filtered/GB119_fixed.vcf.gz
```

done in 3.5 min, with 1368179 records retained (eg 29k removed!) 

one more filter go:

```bash
time ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119_fixed.vcf.gz \
--processes 16 \
--mode phase_change \
--log data/phase_changes/parental.log \
--quality 30 \
--out data/phase_changes/parental/GB119.sam
```

2000 reads in a matter of seconds - it's going to be a long fucking night

```python
# after head -n 300'ing the above file

import classification as rc
import pysam

bam_fname = '../../rcmb/data/phase_changes/parental/GB119_filt.sam'
vcf_fname = '../../rcmb/data/genotyping/vcf_filtered/GB119_fixed.vcf.gz'

pairs = [pair for pair in rc.pairs_creation(bam_fname, vcf_fname)]
```

alright for one - I need to filter out any call with an `*` allele since those are causing
a lot of these

that said I'm still finding some that look like clear good quality matches?? let's up the quality
filter a bit in a few - first I want to look at some bams offline

one finding - based on the read around chromosome 6 at position 7367330 - a
microsat alignment could also give the false appearance of a phase change -
although this one just seems particularly unlucky

for now - now that vcfprep has been updated again to remove the `*` calls, let's
see how many fewer 'phase changes' we get this time

```bash
zcat data/genotyping/vcf_filtered/GB119_fixed.vcf.gz | wc -l # 1368179

time ../readcomb/readcomb/vcfprep.py --vcf data/genotyping/vcf/GB119x3062.vcf.gz \
--no_hets --min_GQ 30 --out data/genotyping/vcf_filtered/GB119_fixed.vcf.gz # 1317318 retained

time ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119_fixed.vcf.gz \
--processes 16 \
--mode phase_change \
--log data/phase_changes/parental.log \
--quality 30 \
--out data/phase_changes/parental/GB119.sam
```

way fewer off the bat - barely up to 850 where I'd normally have 2000 plus

WAIT - after looking at about 15-20 of these on IGV, they seem to consistently be alignment
errors (e.g. around microsats) or occasional sequencing errors. I was thinking about how to
maybe identify certain regions as being problematic for alignment, but turns out - I can just
use these collections of false positives as exactly that! if I get a 'recombination event' 
in the GB119x3062 recombinants that looks exactly the same as one that was 'detected' in GB119
alone, it's likely this was caused by some misalignment - which is a heuristic that saves me countless
hours of reviewing things manually

I think I need to write a workflow for this - but first, pushing the readcomb changes to 0.1.2

next up - regenerating all the VCFs in `data/genotyping/vcf_filtered`

```bash
# from genotyping log
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 16
# reruns readcomb_vcfprep and vcf_tabix
```

once this is done, I'll have to create a workflow for the parental reads - there already
exists `data/phase_changes/parental` for this at least, although it occurs to me
that for each parent A, I'll need to generate a separate bam for each cross it's been
involved with for comparison purposes - which will be a lot of files! 

this was difficult but let's give it a go:

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
```

## 3/1/2022

ran overnight and looks like it just did bamprep for 2.5 hours, but that seems
to have worked without any issues

trying to run filter next:

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16 readcomb_filter
```

things to do once these are done:
- plot read midpoints across the genome - this should give away which regions have more reads than others
- are peaks consistent across samples? this could point to generally messier areas of the genome that
have higher probabilities of false positives in the actual data
- run filter on actual samples! 

## 4/1/2022

so looks like most of the files are just taking 40 min a piece - which means parental filtering
might be done by the end of tomorrow or so

for now, I've got enough 2344 samples that I could start doing some midpoint analyses for these
'phase changes' to see if the same regions of the genome are leading to false positives across
different MT- strains

after this I should also look into why readcomb isn't correctly generating progress bars - has
to do with the `pair_count` attribute - as well as fix up some of the bamprep stderr writing

going to look at midpoints in `false_positives.ipynb`

## 5/1/2022

alright - the parental lines should be done overnight tonight, so I'm going
to queue up the actual recombinants

before I do - need to fix up the progress bar bug in `readcomb`

```bash
time ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--processes 16 --quality 30 --out data/phase_changes/sam/2343x1691.filtered.sam
```

fixed it! was an errant `bai` file check - going to push changes, recompile readcomb, and push to pip

out of curiosity - I'm also going to run a parental sample with quality bumped up to 40 and see if
that leads to any false positives

```bash
time readcomb-filter \
--bam data/alignments/bam_filtered/GB119.sorted.bam \
--vcf data/genotyping/vcf_filtered/GB119x3062.vcf.gz \
--processes 12 \
--mode phase_change \
--quality 40 \
--out GB119_test.sam
```

still some hundreds with GQ40, and they disappear entirely at 50 - I'm still
going to stick to 30, but this is good to know

next up - I'm going to rerun the actual samples, since I also need to do the
mtDNA leakage quantification from those directly - also going to add a rule
that converts the files to BAM at the end of it (why does something tell me I
want to keep them in SAM format? I can't think of any good reasons why right
now)

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 16
```

## 7/1/2022

fixed some weirdness in the CC1691 parental run (see `alignment` log) - here we go again:

```bash
# rerun bamprep on 1691 and start on phase change detection again
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk readcomb_filter --cores 16
```

apparently directly murdering that one garbage corrupted read wasn't enough! time to do it
again I guess????? 

update - all done now (see alignment log - what a nightmare) - rerunning: 

```bash
# run 2: electric boogaloo
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk readcomb_filter --cores 16
```

## 15/1/2022

redoing run on samples now that markers have been updated (although these may be updated once more
afterwards based on the het search) 

```bash
# redoes readcomb_filter and bam_convert for each file
# snakemake actually detected that input files changed this time!
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 16
```

going to create a sorted sam to view in IGV for one of the completed ones - let's do
2932x3062

```bash
time samtools sort -O sam data/phase_changes/bam/2932x3062.filtered.bam \
-o data/phase_changes/sam/2932x3062.sorted.sam
```

going to download this and the corresponding VCF for use with IGV 

## 19/1/2022

today - things to add to readcomb:

- correct indexing of reads w/ indels for query quality comparison (currently just includes deletions)
- code that compares `detection_1` and `detection_2` between the two reads - if a variant is represented more
than once and the calls disagree, that variant should be ignored

alright - changes implemented (among some others) - testing this out before
updating readcomb proper:

```bash
mkdir -p temp_2932 # in root

# the only updates to filter have been the corrected indel quality checks
time python ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/2932x3062.sorted.bam \
--vcf data/genotyping/vcf_filtered/2932x3062.vcf.gz \
--processes 16 \
--out temp_2932/2932x3062.filtered.sam \
--quality 30
```

and then on the parents for false positive detection:

```bash
# half processes since I'm running both overnight concurrently
time python ../readcomb/readcomb/filter.py \
--bam data/alignments/parental_bam_filtered/CC2932.sorted.bam \
--vcf data/genotyping/vcf_filtered/2932x3062.vcf.gz \
--processes 8 --quality 30 \
--out temp_2932/2932x3062.plus.filtered.sam

time python ../readcomb/readcomb/filter.py \
--bam data/alignments/parental_bam_filtered/CC3062.sorted.bam \
--vcf data/genotyping/vcf_filtered/2932x3062.vcf.gz \
--processes 8 --quality 30 \
--out temp_2932/2932x3062.minus.filtered.sam
```

## 21/1/2022

need to sort the three files above to view in IGV:

```bash
# in temp_2932
samtools sort -o 2932x3062.filtered.sorted.bam 2932x3062.filtered.sam
samtools sort -o 2932x3062.plus.filtered.sorted.bam 2932x3062.plus.filtered.sam
samtools sort -o 2932x3062.minus.filtered.sorted.bam 2932x3062.minus.filtered.sam
```

## 23/1/2022

creating a filtered bam that only contains reads with detected phase changes (eg removing
the `no_phase_change` reads) - if this works as expected, will push the readcomb changes

```python
# in readcomb folder - hence weird paths
import pysam
import classification as rc
from tqdm import tqdm

bam_fname = '../../rcmb/temp_2932/2932x3062.filtered.sam'
vcf_fname = '../../rcmb/data/genotyping/vcf_filtered/2932x3062.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)
writer_template = pysam.AlignmentFile(bam_fname, 'r')
writer = pysam.AlignmentFile(
    '../../rcmb/temp_2932/2932x3062.het.filtered.sam', 'wh', template=writer_template) 

for pair in tqdm(reader):
    pair.classify(masking=0, quality=30) # masking size doesn't really matter here
    if 'no_phase_change' not in [pair.call, pair.masked_call]:
        writer.write(pair.rec_1)
        writer.write(pair.rec_2)
# setting masking size to 0 for now since 2932 has some 50 bp reads
```

## 24/1/2022

took several hours of debugging that culminated in a new classification function (`quality_cigar`)
that realigns the query qualities as well - but looks like we're in business, finally! 

took 20 minutes for readcomb to go through the 272k pairs in the filtered sam, and it
only wrote 92k of them - will have to review these in IGV tomorrow

for now, pushing the readcomb changes as 0.2

## 25/1/2022

looking at the output reads in IGV - but first, need to sort:

```bash
# in temp_2932
samtools sort -O bam -o 2932x3062.het.filtered.sorted.sam 2932x3062.het.filtered.sam
samtools index 2932x3062.het.filtered.sorted.bam
```

question - how much of a given chromosome is encompassed by 'false positive' reads?

```python
# in temp_2932
import pysam
from tqdm import tqdm

reader = pysam.AlignmentFile('2932x3062.plus.filtered.sorted.bam', 'rb')

lookup = '0' * 8225626 # length of chr1

for record in tqdm(reader):
    start = record.reference_start
    end = start + sum([basecount for op, basecount in record.cigartuples if op in [0, 2]])
    to_add = '1' * (end - start)
    lookup = lookup[:start] + to_add + lookup[end:]

len(lookup) # 293719
lookup.count('1') # 7931907
# 3.7 percent! not as much as I thought
```

going to get this going across all chromosomes:

```python
# copied from a VCF header and fixed with vim macros lol
lengths = {'chromosome_01': 8225636,
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
'chromosome_17': 6954842}

frac_dict = {}
for chrom in lengths:
    print(chrom)
    lookup = '0' * lengths[chrom]
    for record in tqdm(reader.fetch(chrom)):
        start = record.reference_start
        end = start + sum([basecount for op, basecount in record.cigartuples if op in [0, 2]])
        to_add = '1' * (end - start)
        lookup = lookup[:start] + to_add + lookup[end:]
    frac_dict[chrom] = [lookup.count('1'), lookup.count('0'), lookup.count('1') / lookup.count('0')]
```

just took 3-4 min - looks like rates are consistently ~3% or lower:

```python
>>> for chrom in frac_dict:
...     print(chrom, frac_dict[chrom])
chromosome_01 [293719, 7931917, 0.037030014308016584]
chromosome_02 [314879, 8341005, 0.037750726681017456]
chromosome_03 [318490, 8968404, 0.03551245015278081]
chromosome_04 [131224, 3998849, 0.0328154426436207]
chromosome_05 [119316, 3562844, 0.03348897678371548]
chromosome_06 [296917, 8616442, 0.034459351087142466]
chromosome_07 [240974, 6251133, 0.03854885186413407]
chromosome_08 [153365, 4373618, 0.03506593397045649]
chromosome_09 [203891, 6603257, 0.03087733825898341]
chromosome_10 [272044, 6528203, 0.041672110992872005]
chromosome_11 [133694, 4345828, 0.030763757792531137]
chromosome_12 [349501, 9603238, 0.03639407874718923]
chromosome_13 [210523, 5070915, 0.04151578166859433]
chromosome_14 [148608, 4068695, 0.03652473336045096]
chromosome_15 [68645, 5801998, 0.01183126915934821]
chromosome_16 [282671, 7759804, 0.03642759533617086]
chromosome_17 [247429, 6707413, 0.03688888696730021]
```

implementing a first pass filter for these regions - going to go with a super conservative
heuristic here, where if a read pair overlaps a '1' region at all we ignore it

going to also create a `lengths.csv` file in the root of the `rcmb` project based off of a VCF header

```python
# in temp_2932
import csv
import pysam
from tqdm import tqdm

with open('../lengths.csv', 'r') as f:
    lengths = {line['chrom']: int(line['length']) for line in csv.DictReader(f)}

plus_reader = pysam.AlignmentFile('2932x3062.plus.filtered.sorted.bam', 'rb')
minus_reader = pysam.AlignmentFile('2932x3062.minus.filtered.sorted.bam', 'rb')

# greedy - holding all contig lookups in memory
lookups = {}

for chrom in lengths:
    lookups[chrom] = '0' * lengths[chrom]

for gen in [plus_reader, minus_reader]:
    for chrom in lengths:
        for record in tqdm(gen.fetch(chrom), desc=chrom):
            start = record.reference_start
            end = start + sum([basecount for op, basecount in record.cigartuples if op in [0, 2]])
            to_add = '1' * (end - start)
            lookups[chrom] = lookups[chrom][:start] + to_add + lookups[chrom][end:]
```

since I'm curious about how much including both parents inflates numbers -

```python
frac_dict = {}
for chrom in lengths:
    l = lookups[chrom]
    frac_dict[chrom] = [l.count('1'), l.count('0'), l.count('1') / l.count('0')]

# seems this 'cleanly doubles' the amount of suspect regions...
chromosome_01 [497634, 7728002, 0.06439361687535795]
chromosome_02 [569100, 8086784, 0.07037408195890975]
chromosome_03 [562991, 8723903, 0.06453430305220037]
chromosome_04 [229038, 3901035, 0.05871211101669172]
chromosome_05 [202061, 3480099, 0.05806185398748714]
chromosome_06 [531459, 8381900, 0.06340555244037747]
chromosome_07 [420829, 6071278, 0.06931473077002898]
chromosome_08 [264834, 4262149, 0.062136260370062146]
chromosome_09 [366404, 6440744, 0.05688845884885348]
chromosome_10 [472551, 6327696, 0.07467978866241362]
chromosome_11 [248594, 4230928, 0.05875637685160324]
chromosome_12 [611264, 9341475, 0.06543549064789018]
chromosome_13 [346139, 4935299, 0.07013536565869667]
chromosome_14 [261505, 3955798, 0.06610676278212387]
chromosome_15 [119077, 5751566, 0.020703404950929885]
chromosome_16 [506960, 7535515, 0.06727609194593867]
chromosome_17 [441455, 6513387, 0.06777656540291556]
contig_18 [154, 146574, 0.0010506638285098311]
contig_19 [137, 104692, 0.001308600466129217]
contig_20 [0, 104799, 0.0]
contig_21 [62, 102930, 0.0006023511124064898]
contig_22 [0, 102447, 0.0]
contig_23 [0, 97645, 0.0]
contig_24 [2288, 84094, 0.027207648583727734]
contig_25 [0, 79732, 0.0]
contig_26 [0, 73703, 0.0]
contig_27 [0, 72308, 0.0]
contig_28 [0, 46567, 0.0]
contig_29 [0, 40261, 0.0]
contig_30 [0, 39771, 0.0]
contig_31 [0, 35441, 0.0]
contig_32 [0, 34618, 0.0]
contig_33 [0, 34215, 0.0]
contig_34 [0, 33186, 0.0]
contig_35 [0, 30690, 0.0]
contig_36 [0, 30442, 0.0]
contig_37 [0, 29958, 0.0]
contig_38 [0, 29289, 0.0]
contig_39 [0, 27844, 0.0]
contig_40 [0, 26488, 0.0]
contig_41 [178, 25922, 0.006866754108479284]
contig_42 [0, 25547, 0.0]
contig_43 [0, 22868, 0.0]
contig_44 [0, 21662, 0.0]
contig_45 [0, 21071, 0.0]
contig_46 [0, 20971, 0.0]
contig_47 [0, 19647, 0.0]
contig_48 [0, 18714, 0.0]
contig_49 [0, 18264, 0.0]
contig_50 [0, 17693, 0.0]
contig_51 [0, 17217, 0.0]
contig_52 [0, 16980, 0.0]
contig_53 [0, 14371, 0.0]
contig_54 [0, 14279, 0.0]
contig_55 [0, 12394, 0.0]
contig_56 [0, 10144, 0.0]
contig_57 [0, 7790, 0.0]
cpDNA [10578, 194957, 0.05425811845689049]
mtDNA [449, 15340, 0.029269882659713167]
MT_plus_R [1849, 385381, 0.004797849400982404]
```

and now for iterating through pairs with the lookups - reproducing code from above
just so everything's in one place

```python
# from above
# in temp_2932
import csv
import pysam
from tqdm import tqdm

with open('../lengths.csv', 'r') as f:
    lengths = {line['chrom']: int(line['length']) for line in csv.DictReader(f)}

plus_reader = pysam.AlignmentFile('2932x3062.plus.filtered.sorted.bam', 'rb')
minus_reader = pysam.AlignmentFile('2932x3062.minus.filtered.sorted.bam', 'rb')

# greedy - holding all contig lookups in memory
lookups = {}

for chrom in lengths:
    lookups[chrom] = '0' * lengths[chrom]

for gen in [plus_reader, minus_reader]:
    for chrom in lengths:
        for record in tqdm(gen.fetch(chrom), desc=chrom):
            start = record.reference_start
            end = start + sum([basecount for op, basecount in record.cigartuples if op in [0, 2]])
            to_add = '1' * (end - start)
            lookups[chrom] = lookups[chrom][:start] + to_add + lookups[chrom][end:]

# note - lookups is 114 Mb! yikes! 
# I wonder if I can make a temp tabix file for this in the future, like in all_callables for the salt project
# newly added code:

import readcomb.classification as rc

bam_fname = '2932x3062.het.filtered.sam'
vcf_fname = '../data/genotyping/vcf_filtered/2932x3062.vcf.gz'

reader_template = pysam.AlignmentFile(bam_fname, 'r')
writer = pysam.AlignmentFile('2932x3062.final.sam', 'wh', template=reader_template)

kept, discarded = 0, 0
for pair in tqdm(rc.pairs_creation(bam_fname, vcf_fname)):
    chrom, pos = pair.location.split(':')
    start, end = [int(n) for n in pos.split('-')]
    if '1' in lookups[chrom][start:end]:
        discarded += 1
        continue
    else:
        kept += 1
        writer.write(pair.rec_1)
        writer.write(pair.rec_2)

# 10100 kept, 82601 discarded
```

and now to view these reads in IGV:

```bash
# in temp_2932
samtools sort -O bam -o 2932x3062.final.sorted.bam 2932x3062.final.sam
samtools index 2932x3062.final.sorted.bam
```

more notes to self

- filter by MAPQ in bamprep - there's at least one read pair (chr1:373063) that is clearly misaligned
and has a MAPQ value of 9 and 15 for the two reads
    - what's an appropriate threshold for this?
    - it does seem most 'good' reads have MAPQ 60
- filter reads that overlap but have 'disagreeing' sections - this makes the read suspect
    - `resolve_overlap` does this for variants, but honestly if a read has this maybe we chuck it altogether
- possibly still need to mask entire regions prone to paralogy - the false positive reads work but aren't 100%
    - region at chr1:373133 is a trainwreck in both parents
- probably filter out reads with any `no_match` instances
- need a focal variant attribute - indel proximity seems to matter

also - why is the read at 3354584 a phase change read? review tomorrow

## 26/1/2022

first - reviewing that one weird read - and then some implementations of the things above

```python
import pysam
import readcomb.classification as rc
from tqdm import tqdm

bam_fname = '2932x3062.het.filtered.sam'
vcf_fname = '../data/genotyping/vcf_filtered/2932x3062.vcf.gz'

weird = 'A00516:220:H2NGGDRXY:1:2267:29984:12524'

for pair in tqdm(reader):
    if pair.rec_1.query_name == weird:
        weird_pair = pair
        break

weird_pair.detection
```

ah - looks like there was just a sneaky section with two adjacent SNPs that I missed - that
explains it 

and now for some implementation! 

let's look at the distribution of mapq scores in a given sample bam:

```python
>>> import pysam
>>> from tqdm import tqdm
>>> reader = pysam.AlignmentFile('data/alignments/bam/2932x3062.bam')
>>> scores = {}
>>> counter = 0
... for record in tqdm(reader):
...     if not record.mapq in scores:
...         scores[record.mapq] = 1
...     elif record.mapq in scores:
...         scores[record.mapq] += 1
...     counter += 1
...     if counter > 1e7:
...         break
9978281it [00:19, 520640.93it/s]

>>> for score in sorted(scores.keys()):
...     print(score, scores[score])
0 1469328
1 16634
2 14391
3 20864
4 15535
5 14343
6 18949
7 15732
8 8867
9 20679
10 7633
11 6725
12 10400
13 8911
14 5930
15 12108
16 8205
17 9054
18 11002
19 11012
20 8923
21 17209
22 11665
23 9309
24 13299
25 17428
26 5437
27 25578
28 6391
29 5497
30 9381
31 7090
32 5256
33 12901
34 5268
35 5405
36 6907
37 7347
38 4732
39 11393
40 109752
41 11239
42 13187
43 14799
44 18103
45 24291
46 29701
47 11852
48 11271
49 11491
50 9629
51 10125
52 15178
53 10869
54 16034
55 13437
56 16390
57 18677
58 19002
59 14798
60 7697458
```

super preliminary, but it does tell me there are a lot of reads past MAPQ 40 and even 50

redoing this for a sample parental bam - let's try CC3062

```pythono
>>> for score in sorted(scores.keys()):
...     print(score, scores[score])
0 1218791
1 10694
2 7872
3 13534
4 8170
5 7779
6 12691
7 10389
8 7234
9 14125
10 7645
11 6177
12 10875
13 8957
14 7264
15 13144
16 8001
17 9276
18 13618
19 13000
20 11381
21 14152
22 15169
23 12206
24 20176
25 25138
26 4029
27 46745
28 4259
29 3829
30 9880
31 5481
32 3987
33 11142
34 4388
35 4105
36 9313
37 6196
38 4061
39 12068
40 123770
41 7193
42 13075
43 8909
44 10658
45 15449
46 11709
47 14588
48 21505
49 20310
50 27796
51 9664
52 11435
53 8135
54 20870
55 11095
56 9464
57 16389
58 10316
59 13781
60 7976949
```

looks like lots of reads above MAPQ 40 and even 50 here as well - let's try 40 as a
preliminary filter then

alright - implemented MAPQ filter into `readcomb-filter` and I've added in new overlap
attributes to `classification`. the final thing to add is a `focal_variants` attribute
and maybe an indel proximity (to focal variants) one after that but _specifically for
gene conversions_ since that's where this seems to be an issue most

## 27/1/2022

done all of the above and they seem to function fine - going to give the mapq
functionality a test drive

```bash
time python ../readcomb/readcomb/filter.py \
--bam data/alignments/bam_filtered/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--min_mapq 40 --processes 16 --quality 30 --out 2343x1691.mapq.test.sam
```

if this looks good in the morning - push readcomb changes to 0.2.1 and
redo `readcomb-filter` across all sample bams and parental bams (for false positives)

## 30/1/2022

readcomb's now at 0.2.1 - time to redo phase change detection with the mapq
filters in place

```bash
# min mapq set to 40
# being run with VCFs that were generated with GQ20 filter
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk \
--cores 16
# will take a full day
```

having a look at a file - first, removing all the no phase change reads
with readcomb:

```python
# in a new root folder called temp_2344
import readcomb.classification as rc
from tqdm import tqdm
import pysam

bam_fname = 'data/phase_changes/bam/2344x1952.filtered.bam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x1952.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)

writer = pysam.AlignmentFile('temp_2344/2344x1952.chr1.bam', 'wh',
    template=pysam.AlignmentFile(bam_fname, 'rb'))

total = 0
found = 0
for pair in tqdm(reader):
    total += 1
    if pair.rec_1.reference_name == 'chromosome_01':
        pair.classify(masking=0)
        if pair.call != 'no_phase_change':
            found += 1
            writer.write(pair.rec_1)
            writer.write(pair.rec_2)

# total 595309, found 40931
```

and then prepping this file for IGV:

```bash
time samtools sort -o temp_2344/2344x1952.chr1.sorted.bam temp_2344/2344x1952.chr1.bam
samtools index temp_2344/2344x1952.chr1.sorted.bam
```

## 1/2/2022

happy 'it's not January anymore' day! 

today - need to implement a new quality metric discussed with Rob that should
catch paralogous reads - basically involves counting the number of mismatches
and dividing by the number of SNPs in the same region (since if there's a large
number of mismatches and a low number of SNPs, that means those sites likely did
not pass SNP filters, which is reflective of a 'crap' region)

other things to work on could include quantifying gene drive from the full bams - 
but that's another thing to look at, in another folder with another log

for now - implementing unexpected matches direct into `readcomb.classification` and
having a look at a few reads based off that 

## 6/3/2022

where the hell did Feb go? (you had medical problems including two weeks of having your
digestive system destroyed by an antibiotic. that's where it went)

looks like I did incorporate said metric - getting the mount going and then looking at
`temp_2344/2344x1952.chr1.sorted.bam` - and also filtering by the new metrics I added

time to get some sample COs to look at - filtering further still from the first pass

```python
import readcomb.classification as rc
import pysam
from tqdm import tqdm
import time

bam_fname = 'temp_2344/2344x1952.chr1.bam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x1952.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)
writer = pysam.AlignmentFile('temp_2344/2344x1952.chr1.COs.sam', 'wh',
    template=pysam.AlignmentFile(bam_fname))

for pair in tqdm(reader):
    pair.classify(masking=0)
    if pair.min_variants_in_haplotype:
        if pair.call == 'cross_over' and pair.min_variants_in_haplotype > 1:
            writer.write(pair.rec_1)
            writer.write(pair.rec_2)
            time.sleep(0.01) # since it sometimes doesn't write correctly
```

sorting and indexing for IGV:

```bash
samtools sort -O bam -o 2344x1952.chr1.COs.bam 2344x1952.chr1.COs.sam
samtools index 2344x1952.chr1.COs.bam
```

## 22/4/2022

alright new parental sequences just dropped, fire emoji 100 emoji

now that variant calling is done, trying out one sample:

```bash
mkdir new_test

time readcomb-filter \
--bam data/alignments/bam_filtered/2344x2931.sorted.bam \
--vcf data/genotyping/vcf_filtered/2344x2931.vcf.gz \
--min_mapq 40 --processes 16 --quality 30 --out new_test/2344x2931.filtered.sam
```

C-c'd out of it around 4 min in, after 1.3m pairs had been processed

```bash
# in new_test - get first 2000 read pairs
head -n 2070 2344x2931.filtered.sam > 2344x2931.sam
rm 2344x2931.filtered.sam
mv -v 2344x2931.sam 2344x2931.filtered.sam

# prep for viewing in IGV
samtools sort -O bam -o 2344x2931.filtered.bam 2344x2931.filtered.sam
samtools index 2344x2931.filtered.bam
```

separately pulling up IGV - but let's get a list of crossovers first:

```python
import readcomb.classification as rc
from tqdm import tqdm
import pysam

bam_fname = 'new_test/2344x2931.filtered.sam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x2931.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)

cos = []
for pair in tqdm(reader):
    pair.classify(masking=0)
    if pair.min_variants_in_haplotype and pair.call == 'cross_over':
        if pair.min_variants_in_haplotype > 1:
            cos.append(pair)
            if len(cos) == 10:
                break
```

to incorporate:
- set numerical values in classification.Pair to -1 by default, not None
- add an indel proximity thing to crossovers as well - proximity to midpoint basically
- most importantly - add a 'proximity to end of read' thing for CO breakpoints 
    - should explain some false positives

but first - going to get the parental false positives workflow going - after that, will also 
need to incorporate some stuff to filter out recombination events that have actually been
caused by het calls

in the meantime - getting the parental phase changes going:

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
```

## 24/4/2022

today:

- look at parental phase change counts, out of curiosity
- up the quality filters (eg min variants in haplotype) and find 'better' calls
- debug the weird div by zero error that's happening in `mismatch_variant_ratio`
- plan out lab meeting presentation 

re parental phase changes - they're currently sorted by read name, but how do I
want to handle these? it may be better to have them sorted by position if I want
to create lookups - but having them by read name lets me ID problematic conversion tracts
and misleading variants

I could also use these to filter VCFs directly - need to talk to Rob about this one

meanwhile - it seems that upping the min variants in haplotype filter to 3 just means
finding more false positives. the common threads seem to be:

- very short haplotype tract for one of the two parents
- two read types in the parents

wait - I forgot to update the snakefile to actually do more than just readcomb-bamprep...

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
# should also readcomb-filter
```

## 27/4/2022

so there's been some movement on changing up the genotyping workflow - need to
finalize switching over to freebayes, but in the meantime, I need to figure out
how to create a lookup for false positives in the parental read sets

basically - if a phase change is called on a parental read, we know that that phase change
is 'a lie' - and we want to filter out any supposed recombinant reads with that same phase change

here's an example:

```python
>>> false_pair.detection
[('2', 3424823, 'C'), ('2', 3424875, 'C'), ('1', 3425013, 'T'), 
('1', 3425020, 'G'), ('1', 3425025, 'T'), ('2', 3425054, 'A'), 
('2', 3425089, 'C'), ('2', 3425105, 'A'), ('2', 3425123, 'A'), 
('2', 3425148, 'C'), ('2', 3425188, 'G'), ('2', 3425205, 'T'), 
('2', 3425223, 'T'), ('2', 3425242, 'A'), ('2', 3425244, 'A')]

>>> cos[2].detection
[('1', 3425013, 'T'), ('1', 3425020, 'G'), ('1', 3425025, 'T'), 
('2', 3425089, 'C'), ('2', 3425105, 'A'), ('2', 3425123, 'A'), 
('2', 3425148, 'C'), ('2', 3425188, 'G'), ('2', 3425205, 'T'), 
('2', 3425223, 'T'), ('2', 3425242, 'A'), ('2', 3425244, 'A'), 
('2', 3425281, 'T'), ('2', 3425328, 'G'), ('2', 3425367, 'G'),
('2', 3425371, 'G'), ('2', 3425410, 'G'), ('2', 3425434, 'G'), 
('2', 3425445, 'A'), ('2', 3425449, 'T')]
```

where `cos[2]` is the detection output for a false phase change, while
`false_pair` is the equivalent parental read that has that false event

thought: could use tabix and create bed intervals with phase changes! 

so something like 

```
left_bound right_bound left_parent right_parent
3424875 3425013 2 1 
3425025 3425054 1 2
```

that could then be checked against each individual read! 

## 3/5/2022

today - implementing this in a new readcomb tool called `false_positives`

should take in -

1. readcomb output for a given cross
2. readcomb output on plus parent bam
3. readcomb output on minus parent bam
4. which filtering method to use (midpoint/overlap/nuclear)
5. bed files of false positives (if appropriate)
6. log file
7. file to write to

## 12/5/2022

day 5 of covid right now but I think I got it done! let's
give it a go - first, would need to run readcomb filter on parents,
and then generate bed files

```bash
python ../readcomb/readcomb/false_positives.py \
--fname freebayes-test/2344x2931.filtered.sam \
--false_plus freebayes-test/2344x2931.plus.FP.sam \ # ran readcomb-filter on parent files for these
--false_minus freebayes-test/2344x2931.minus.FP.sam \
--vcf freebayes-test/2344x2931.prepped.vcf.gz \
--method midpoint \
--false_bed_out 2344x2931.test.bed \
--out freebayes-test/2344x2931.fpd.sam
```

## 14/5/2022

alright, this script has now been debugged after several hours and runs just fine on the 
three filter modes - also logs things now too

time to regen the parental false positive sams with the new VCFs and then it's off to the races

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
# regens filtered sams in data/phase_changes/parental/
```

## 15/5/2022

done in a day - now to get a workflow going to filter out the false positives 

I should first look at some of the filtered reads returned in `freebayes-test` though
to see which of the three methods I want to use

## 16/5/2022

some notes from yesterday:

- there's still a bug in midpoint calculation that I've since fixed (the
start parameter in relative midpoint calc was incorrect)
- midpoint calculation is wrong for gene conversions - needs to be between the focal variants 

now that both of these are fixed, going to create a giant data frame of calls to look at

```python
import csv
import readcomb.classification as rc
from tqdm import tqdm

# picking midpoint bc there was one false GC filtered in midpoint but not overlap
bam_fname = 'freebayes-test/fpd/2344x2931.fpd.midpoint.sam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x2931.vcf.gz'
reader = rc.pairs_creation(bam_fname, vcf_fname)

fieldnames = [
    'chromosome', 'midpoint', 'rel_midpoint', 'call', 'masked_call', 'mask_size', 'var_count', 
    'min_vars_in_hap', 'var_skew', 'mismatch_var_ratio', 'var_per_hap', 'gc_length', 
    'indel_proximity', 'detection', 'read_name']

with open('freebayes-test/2344x2931.midpoint.tsv', 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for pair in tqdm(reader):
        pair.classify(masking=20)
        writer.writerow({
            'chromosome': pair.rec_1.reference_name, 'midpoint': pair.midpoint,
            'rel_midpoint': pair.relative_midpoint, 'call': pair.call, 'masked_call': pair.masked_call,
            'mask_size': 20, 'var_count': len(pair.variants_filt),
            'min_vars_in_hap': pair.min_variants_in_haplotype,
            'var_skew': pair.variant_skew, 'mismatch_var_ratio': pair.mismatch_variant_ratio,
            'var_per_hap': pair.variants_per_haplotype, 'gc_length': pair.gene_conversion_len,
            'indel_proximity': pair.indel_proximity,
            'detection': pair.detection, 'read_name': pair.rec_1.query_name})
```
     
## 17/5/2022

today - running readcomb-filter on all crosses and then employing midpoint filtering on all

updating `phase_change_detection.smk` to do both `readcomb-filter` and `readcomb-fp`
in one go - will also generate sorted bams separately for IGV viewing as necessary

here goes:

```bash
# readcomb-filter, readcomb-fp, and bam creation
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 16
```

## 19/5/2022

now that readcomb filter has output a few files, going to make equivalent data
frames for some of these other crosses and look at crossovers + gene conversions:

making a dir in `data/phase_changes` called `event_summaries` - I should
maybe update the phase change snakemake workflow with this at the tail end
at some point - perhaps with a script attached

```python
import csv
import readcomb.classification as rc
from tqdm import tqdm

# picking midpoint bc there was one false GC filtered in midpoint but not overlap
bam_fname = 'data/phase_changes/sam/2343x1952.filtered.sam'
vcf_fname = 'data/genotyping/vcf_filtered/2343x1952.vcf.gz'
reader = rc.pairs_creation(bam_fname, vcf_fname)

fieldnames = [
    'chromosome', 'midpoint', 'rel_midpoint', 'call', 'masked_call', 'mask_size', 'var_count', 
    'min_vars_in_hap', 'var_skew', 'mismatch_var_ratio', 'var_per_hap', 'gc_length', 
    'indel_proximity', 'detection', 'read_name']

with open('data/phase_changes/event_summaries/2343x1952.midpoint.tsv', 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for pair in tqdm(reader):
        pair.classify(masking=50)
        writer.writerow({
            'chromosome': pair.rec_1.reference_name, 'midpoint': pair.midpoint,
            'rel_midpoint': pair.relative_midpoint, 'call': pair.call, 'masked_call': pair.masked_call,
            'mask_size': 50, 'var_count': len(pair.variants_filt),
            'min_vars_in_hap': pair.min_variants_in_haplotype,
            'var_skew': pair.variant_skew, 'mismatch_var_ratio': pair.mismatch_variant_ratio,
            'var_per_hap': pair.variants_per_haplotype, 'gc_length': pair.gene_conversion_len,
            'indel_proximity': pair.indel_proximity,
            'detection': pair.detection, 'read_name': pair.rec_1.query_name})
```

## 24/5/2022

today - continuing to diagnose what's weird with 3071 and 3086 (see alignment
and genotyping logs) 

I've regen'd the VCFs using older 2 x 150 bp sequences - now to rerun the
parental false positive file workflow and then the phase change workflow 

```bash
mkdir data/phase_changes/parental_weird
mv -v data/phase_changes/parental/3071* data/phase_changes/parental_weird
mv -v data/phase_changes/parental/3086* data/phase_changes/parental_weird

time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 4
```

## 26/5/2022

looks like the 3071 file switch dramatically dropped the number of recombination
events to something that looks a lot more normal - which is good

let's create a set of quality crossovers for each - going
to create a simple filtering script since I anticipate re-running this often

actually - I should do this in R instead - that way I should be able to plot and
query things more easily as well

## 29/5/2022

alright, redid some variant calls (see genotyping log) - now to remake
the false positives and get phase changes for each sample

```bash
rm -v data/phase_changes/parental/*sam
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
# should take about ~22 hours
```

## 30/5/2022

and now for the phase change runs:

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 16
```

## 1/6/2022

the more things change...

the VCFs have been re-regenerated - now for `readcomb-fp` and the phase change runs
(again)

```bash
rm -v data/phase_changes/parental/*sam
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 16
```

## 3/6/2022

and now for phase changes... again

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 20
```

this is looking quite slow so far, likely since there are so many more new variants...

took 34 hours in total - would have taken 21 days on a single core! 
back to RStudio to get some sample reads out

## 7/6/2022

one more update to summarise cross - going to add read lengths as another column
to filter off

I'll likely have to redo the summarising after updating the masked call algorithm,
but this should be handy for the time being

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 20
```

## 9/6/2022

going to separately do 'nuclear' filtering as well to see if that makes a difference,
in the quest to reduce min vars in hap to 1 and still have good quality data via other
filters

modifying the last three rules of `phase_change_filter` temporarily for this - will
output files in a new subdir of `phase_changes` called `nuclear`

```bash
mkdir -p data/phase_changes/nuclear

# server is clear right now - upping threads but lowering priority
nice -n 5 time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 24
# took 34 hours
```

## 9/7/2022

getting back on this post evol2022 + salt stuff

need to get expected CO/NCO counts - for that I need:

1. the number of kept (filtered) reads per cross
2. the 'effective sequence space' across these reads

this second thing is to account for the fact that a LOT
of these reads seem to overlap, which is an issue

I should be able to get both these things by
just getting a distribution of sequence space counts - eg
something like

```
effective_length count
1 100
2 0
3 0
4 12
5 15
```

all the way up to ~550 - this way, the sum of the second
column simply gives me the number of kept reads

don't think I need a script for this necessarily - let's just try this in a console:

```python
import csv
from tqdm import tqdm
import readcomb.classification as rc

with open('data/alignments/samples.txt', 'r') as f:
    samples = [s.rstrip('\n') for s in f]

with open('data/alignments/read_lengths.tsv', 'a', newline='') as f:
    fieldnames = ['sample', 'effective_length', 'count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for sample in tqdm(samples, desc='samples'):
        bam_fname = f'data/alignments/bam_prepped/{sample}.sorted.bam'
        vcf_fname = f'data/genotyping/vcf_filtered/{sample}.vcf.gz'

        counts = {k: 0 for k in range(1000)}

        reader = rc.pairs_creation(bam_fname, vcf_fname)
        for pair in tqdm(reader, desc=f'{sample}'):
            start = pair.rec_1.reference_start
            end = pair.rec_2.reference_start + len(pair.segment_2)
            effective_length = end - start
            try:
                counts[effective_length] += 1
            except:
                counts[effective_length] = 1
        
        for i in sorted(counts.keys()):
            writer.writerow({
                'sample': sample,
                'effective_length': i,
                'count': counts[i]})

```

I anticipate this will take about 6.5 min per sample from the test
run - which makes for 3.5 hours for all crosses

after that's done, will first need to transform (possibly using
R) the data into a per sample total count to get the number of reads,
and then will get the total number of 'effective sequence space' we have
in bp to work with

## 10/7/2022

today - getting effective sequence space for each cross:

```R
library(tidyverse)

d = read_tsv('data/alignments/read_lengths.tsv', col_types = cols())

eff_seq = d %>%
    mutate(effective_sequence = effective_length * count) %>%
    group_by(sample) %>%
    summarise(sum_seq = sum(effective_sequence))

write_tsv(eff_seq, 'data/alignments/effective_sequence.tsv')

read_counts = d %>%
    group_by(sample) %>%
    summarise(read_count = sum(count))

write_tsv(read_counts, 'data/alignments/read_counts.tsv')
```

also going to get variant counts per cross while I'm at it -
I could probably do this with a bash oneliner involving `wc -l`
but my brain isn't working enough for that right now

```python
import csv
from cyvcf2 import VCF
from tqdm import tqdm

with open('data/alignments/samples.txt', 'r') as f:
    samples = [s.rstrip('\n') for s in f]

with open('data/genotyping/variant_counts.tsv', 'w', newline='') as f:
    fieldnames = ['sample', 'var_count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    for sample in tqdm(samples, desc='samples'):
        reader = VCF(f'data/genotyping/vcf_filtered/{sample}.vcf.gz')
        counter = 0
        for record in tqdm(reader, desc=f'{sample}'):
            counter += 1
        writer.writerow({
            'sample': sample,
            'var_count': counter
            })

```

## 13/7/2022

readcomb debugging notes after going through that - 

- make `get_midpoint` return both the midpoint and relative midpoint in line 763
- debug `overlap_disagree` 
- reduce base quality to 20 - this could happen in `summarise_cross`
- add MAPQ filter (<50) - this could be independent of `rc.classification`, in `summarise_cross`
- fix no match bug (reads with one hap + no match called crossover)
- incorporate 'bounds of read' metric (proximity to read end)

all of these were pretty straightforward, save for the last one... going to 
update `Pair._describe()` 

recompiling readcomb now that these are done - and off we go again tomorrow! 

## 14/7/2022

today - updating `summarise_cross` as specified above and rerunning with new readcomb

- adding 'min end proximity' as a column
- adding MAPQ > 50 filter
- reducing classify base quality to 20

redoing `summarise_cross` - 

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 20
```

## 2/8/2022

today - seeing whether implementing the above filters on the midpoint method
set will do the trick

updating `phase_change_detection.smk` based on changes to `phase_change_nuclear`
(eg new metrics, mapq filter) and rerunning:

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 15
```

looks like 2344x1952 was bonked - regenerating this while I look at the number of crossovers
across the other crosses:

```R
# in event_summaries
library(tidyverse)
library(fs)

fnames = dir_ls('.', regexp = '.*75\\..*')

d_all = map(
  fnames,
  ~ read_tsv(., col_types = cols()) %>% 
    filter(masked_call == 'cross_over') %>% 
    mutate(
      indel_proximity = as.numeric(
        ifelse(indel_proximity == 'N/A', Inf, indel_proximity)
      )
    )
)

names(d_all) = str_extract(fnames, '[GB0-9]{4,5}x[0-9]{4}')

rcmb_qual = d_all %>% 
  map_dfr(
    ~ filter(.,
      call == 'cross_over',
      min_vars_in_hap >= 1, mismatch_var_ratio <= 1.5, 
      outer_bound >= 0.1, outer_bound <= 0.9,
      read1_length >= 150, read2_length >= 150, effective_length >= 200,
      indel_proximity >= 5) %>% 
      nrow()
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = 'cross',
    values_to = 'count') %>% 
  arrange(desc(count)) %>% 
  mutate(s = sum(count))
```

looks like marginal increases (20-40 COs) over the nuclear set... I'm curious whether
those 20-40 COs are actually legit too, given that they must have been removed from the
nuclear set since false positives were nearby


## 4/8/2022

today - trying to check on nuclear filtering, and make sure it's doing its job right

after eyeballing these bed files and parental fp files, these do look right - but
it occurs to me that these are all reads that just passed `filter`, and many might
still be `no_phase_change` or `complex` calls

prepped one of the fp sets as a read name sorted sam for use with `pairs_creation`:

```python
>>> import readcomb.classification as rc
>>> from tqdm import tqdm
>>> event_counts = {}
>>> bam_fname = '2344x1952.fp.sort.sam'
>>> vcf_fname = '../genotyping/vcf_filtered/2344x1952.vcf.gz'
>>> reader = rc.pairs_creation(bam_fname, vcf_fname)
>>> for pair in tqdm(reader):
...     pair.classify(masking=0)
...     if pair.call not in event_counts:
...         event_counts[pair.call] = 1
...     elif pair.call in event_counts:
...         event_counts[pair.call] += 1
>>> event_counts
{'no_phase_change': 4536, 'gene_conversion': 12550, 'complex': 2449, 'cross_over': 13588}
>>> (event_counts['no_phase_change'] + event_counts['complex']) / sum(event_counts.values())
0.210880656945325
```

20% of reads! 

let's try removing these, writing to file, and then seeing how many crossovers are left

```python
import pysam

bam_reader = pysam.AlignmentFile('../alignments/parental_bam/CC2344.bam', 'r')
bam_writer = pysam.AlignmentFile('2344x1952.plus.cleaned.sam', 'w', template=bam_reader)

reader = rc.pairs_creation(bam_fname, vcf_fname)
counter = 0
kept = 0
for pair in tqdm(reader):
    counter += 1
    pair.classify(masking=0)
    if pair.call in ['cross_over', 'gene_conversion']:
        kept += 1
        bam_writer.write(pair.rec_1)
        bam_writer.write(pair.rec_2)
    else:
        continue
# 33125 total, 26139 kept
```

repeating with minus:

```bash
# in data/phase_changes - setting up for the above
samtools view -O sam parental_bam/2344x1952.minus.filtered.sam.bam > 2344x1952.minus.sam
cat header 2344x1952.minus.sam > 2344x1952.minus.h.sam
samtools sort -n -O sam 2344x1952.minus.h.sam > 2344x1952.minus.sam
```

and then back to python:

```python
bam_reader = pysam.AlignmentFile('../alignments/parental_bam/CC1952.bam', 'r')
bam_writer = pysam.AlignmentFile('2344x1952.minus.cleaned.sam', 'w', template=bam_reader)

bam_fname = '2344x1952.minus.sam'
reader = rc.pairs_creation(bam_fname, vcf_fname)
counter = 0
kept = 0
for pair in tqdm(reader):
    pair.classify(masking=0)
    if pair.call in ['cross_over', 'gene_conversion']:
        kept += 1
        bam_writer.write(pair.rec_1)
        bam_writer.write(pair.rec_2)
    else:
        continue
# 18135 total, 14169 kept
```

and now let's generate a bed file from these - going to move them to a folder
called `nuclear_test`

```bash
# in project root

# going to run for 2-3 hours
time readcomb-filter --bam data/alignments/bam_prepped/2344x1952.sorted.bam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--processes 20 --log data/phase_changes/phase_change_filter.log \
--quality 30 --min_mapq 40 --out data/phase_changes/nuclear_test/2344x1952.init.sam

# need to sort the two filtered fp sams by position, compress them, and index
# running these in nuclear_test
samtools sort -O bam 2344x1952.plus.cleaned.sam > 2344x1952.plus.cleaned.bam
samtools index 2344x1952.plus.cleaned.bam
samtools sort -O bam 2344x1952.minus.cleaned.sam > 2344x1952.minus.cleaned.bam
samtools index 2344x1952.minus.cleaned.bam

time readcomb-fp --fname data/phase_changes/nuclear_test/2344x1952.init.sam \
--false_plus data/phase_changes/nuclear_test/2344x1952.plus.cleaned.bam \
--false_plus data/phase_changes/nuclear_test/2344x1952.minus.cleaned.bam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--method nuclear --false_bed_out data/phase_changes/nuclear_test/2344x1952.nuclear.bed.gz \
--log data/phase_changes/false_positives.log \
--out data/phase_changes/nuclear_test/2344x1952.filtered.sam
```

while this runs - I'm curious as to how much sequence space is blocked off
for each cross by the nuclear filtering, and how much variation there is on that front

```python
# in data/phase_changes/nuclear/fp_bed
import os
import csv
import gzip
from glob import glob
from tqdm import tqdm

fnames = glob('*bed.gz')

seq_space_count = {}

for fname in tqdm(fnames, desc='crosses'):
    cross = fname.rstrip('.nuclear.bed.gz')
    seq_space_count[cross] = {}
    with gzip.open(fname, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for line in tqdm(reader, desc=cross):
            chrom, start, end = line['#chrom'], int(line['start']), int(line['end'])
            if chrom not in seq_space_count[cross]:
                seq_space_count[cross][chrom] = end - start
            else:
                seq_space_count[cross][chrom] += end - start

for cross in seq_space_count:
    print(cross, sum(seq_space_count[cross].values()))

# I should save this

with open('nuclear_blocked.tsv', 'w') as f:
    fieldnames = ['cross', 'blocked_seq']
    fieldnames.extend([f'chromosome_0{i}' for i in range(1, 10)])
    fieldnames.extend([f'chromosome_{i}' for i in range(10, 18)])
    fieldnames.extend(['mtDNA', 'cpDNA'])
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    for cross in seq_space_count:
        out_dict = {'cross': cross, 'blocked_seq': sum(seq_space_count[cross].values())}
        for scaffold in fieldnames[2:]:
            try:
                out_dict[scaffold] = seq_space_count[cross][scaffold]
            except KeyError:
                out_dict[scaffold] = 0
        writer.writerow(out_dict)

```
                
## 5/8/2022

today - how do the above numbers change if just COs/GCs are kept in the nuclear filters?

but first, running this (originally written up yesterday) to see how the actual
CO counts and such are impacted:

```bash
time readcomb-fp --fname data/phase_changes/nuclear_test/2344x1952.init.sam \
--false_plus data/phase_changes/nuclear_test/2344x1952.plus.cleaned.bam \
--false_minus data/phase_changes/nuclear_test/2344x1952.minus.cleaned.bam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--method nuclear --false_bed_out data/phase_changes/nuclear_test/2344x1952.nuclear.bed.gz \
--log data/phase_changes/false_positives.log \
--out data/phase_changes/nuclear_test/2344x1952.filtered.sam
# took 14 min
```

once the bed file is made, I'll do a similar calculation to the thing above just for 2344x1952
to then compare against the values in `nuclear_blocked.tsv`

here goes:

```python
import os
import csv
import gzip
from tqdm import tqdm

chroms = {}

with gzip.open('data/phase_changes/nuclear_test/2344x1952.nuclear.bed.gz, 'rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in tqdm(reader, desc=cross):
        chrom, start, end = line['#chrom'], int(line['start']), int(line['end'])
        if chrom not in chroms:
            chroms[chrom] = end - start
        else:
            chroms[chrom] += end - start

>>> sum(chroms.values())
7809839

# from nuclear_blocked.tsv
>>> 7809839 / 10876012
0.71807929229942
```

next up - does the amount of crossovers scale with how much sequence is blocked off? checked
this in RStudio after exporting but it doesn't seem that's the case

how many COs do I get off of the new sam filtered with a smaller fp set? 

```python
import readcomb.classification as rc
from tqdm import tqdm

bam_fname = 'data/phase_changes/nuclear_test/2344x1952.filtered.sam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x1952.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)

counter = 0
for pair in tqdm(reader):
    pair.classify(masking=75, quality=20)
    
    if pair.call == 'no_phase_change':
        continue
    if all([pair.rec_1.mapq >= 50, pair.rec_2.mapq >= 50]) and pair.call == 'cross_over':
        counter += 1

# 19537! 
```

but looks like `summarise_cross` mistakenly filters things of mapq 50 - which is probably
a lot! let's see how many are in the original with those reads included:

```python
# rerunning the above with this incorrect filter
counter = 0
for pair in tqdm(reader):
    pair.classify(masking=75, quality=20)
    
    if pair.call == 'no_phase_change':
        continue
    if all([pair.rec_1.mapq > 50, pair.rec_2.mapq > 50]) and pair.call == 'cross_over':
        counter += 1
# 19475

# and now trying it with the original
bam_fname = 'data/phase_changes/nuclear/sam/2344x1952.filtered.sam'
reader = rc.pairs_creation(bam_fname, vcf_fname)

old_counter = 0
for pair in tqdm(reader):
    pair.classify(masking=75, quality=20)
    
    if pair.call == 'no_phase_change':
        continue
    if all([pair.rec_1.mapq >= 50, pair.rec_2.mapq >= 50]) and pair.call == 'cross_over':
        old_counter += 1
# 14234 - why so much more than in the existing results?
```

going to rerun `summarise_cross.py` at the command line for this cross after changing
the mapq thing and removing the requirement for `no_phase_change` to not be in `masked_call`

```bash
time python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/nuclear/sam/2344x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--mask_size 75 --remove_uninformative --out mapq_filter_test.tsv

# and then the with the cleaned fp
time python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/nuclear_test/2344x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--mask_size 75 --remove_uninformative --out mapq_filter_test.tsv

# crossover calls went from 891 to 1386!
```

this jump alone makes it worth it to apply something similar to readcomb-fp I think...
and it makes sense - I should only be filtering reads overlapping with false positives
that are actually called as crossovers or gene conversions, instead of anything indiscriminately
passed by `readcomb-filter` (since many reads end up getting 'failed' on base quality etc after)

couple things now - I need to update classification so that pairs without masked calls
still have quality metrics calculated for them (since there's a lot more of those than expected)
and then rerun `summarise_cross`

I should also add these things to `summarise_cross`/readcomb to give the
logistic regression more predictors:

- read mapq (summarise cross)
- some measure of (focal?) variant quality (min variant GQ? mean? idk) - do in summarise cross
- maybe even use the QUAL scores from the VCF directly
- number of indels in the read (eg insertions/deletions in the alignment) - this is a readcomb addition
- size of most proximate indel, if there is one

and fix the line in `Pair._describe()` where quality metrics aren't calculated
if the masked call isn't present and/or is `no_phase_change`

## 7/8/2022

continuing work on this - didn't get to finish it earlier

I've added `proximate_indel_length` and `indels` attributes to readcomb,
but still need to modify `summarise_cross` with read mapq and some kind of variant
quality metric (e.g. based off of `pair.variants_filt`)

going to test the readcomb changes before recompiling:

```python
import classification as rc
from tqdm import tqdm

bam_fname = '../../rcmb/data/phase_changes/nuclear/sam/2344x1952.filtered.sam'
vcf_fname = '../../rcmb/data/genotyping/vcf_filtered/2344x1952.vcf.gz'
reader = rc.pairs_creation(bam_fname, vcf_fname)
pairs = []
for pair in tqdm(reader):
    pair.classify(masking=75)
    if pair.indels:
        pairs.append(pair)
    if len(pairs) == 10:
        break

```

looks good - going to recompile and then update `summarise_cross` as needed

added the following metrics:

- `min_var_qual` - lowest QUAL value in variants spanned by read
- `min_var_depth` - lowest allele depth in variants spanned by read
- `avg_diff_parental_gq` - average of pairwise differences in GQ of parental alleles
- `avg_hap_var_gq` - average GQ of alleles matched by recombinant haplotype
- `avg_hap_var_depth` - average GQ of alleles matched by recombinant haplotype

these metrics need to disregard no match variants - so those are now removed
from consideration just for calculating these metrics

and now almost ready to run this again on the nuclear data! before I do, let's sum up the changes:

- added new metrics (see above) both in readcomb and in `summarise_cross` directly
- fixed bug in readcomb that should now mean metrics are still calculated if masked call is `no_phase_change`
- mapq filtering in `summarise_cross` is now >= 50, not just > 50

the main outstanding thing left to do is fix up `readcomb-fp` to only keep fp events
that are actually crossovers or noncrossovers, but before that let's just make sure
`summarise_cross` is working as intended:

```bash
python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/nuclear/sam/2344x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--mask_size 75 --mapq 50 --remove_uninformative --out sc_test.tsv
```

this looks good, though I forgot to fix the issue in midpoint assignment where
it relies on the masked call still.. one more reason to recompile I suppose

also I just realized that all the 'cleaned' stuff I did above was wrong cause
I was using the original parental bams and not parental bams that had been
bamprepped... brb about to swan dive into the ocean

need to test how many of the retained false plus/minus reads are actually
`no_phase_change` bc readcomb-fp should be filtering these already! 

still - it looks like some still slipped through, likely due to prev bugs
in readcomb that existed when I first did the fp runs:

```python
>>> import readcomb.classification as rc
>>> from tqdm import tqdm
>>> bam_fname = 'test.sorted.sam' # 2344x1952 plus fp reads, read name sorted
>>> vcf_fname = '../../genotyping/vcf_filtered/2344x1952.vcf.gz'
>>> reader = rc.pairs_creation(bam_fname, vcf_fname)
>>> d = {}
>>> total = 0
... for pair in tqdm(reader):
...     total += 1
...     pair.classify(masking=0)
...     if pair.call in d:
...         d[pair.call] += 1
...     elif pair.call not in d:
...         d[pair.call] = 1
33125it [04:03, 136.27it/s]
>>> total
33125

>>> d
{'no_phase_change': 4536, 'gene_conversion': 12550, 'complex': 2450, 'cross_over': 13589}
```

going to redo this and also set some read length filters - the reads should
both be at least 100 bp long and specifically not be `no_phase_change` calls

going to recompile readcomb as 0.4.0, and then will first run `summarise_cross`
on existing filtered files just to have something for tomorrow's meeting in the meantime
(regenerating all the fp files will take a while)

```bash
# after clearing `event_summaries` - this should just run the final summarise_cross rule
time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 4
```

and now to regen basically everything overnight:

```bash
rm -v data/phase_changes/nuclear/fp_bed/*bed.gz*tbi
rm -v data/phase_changes/nuclear/sam/*.filtered.sam
rm -v data/phase_changes/nuclear/bam/*filtered*bam
time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 30
```

## 9/8/2022

today - lots of things...

first - let's dial back on basically all filters we have going here

that includes the mapq filter in `readcomb-filter` and the base quality filters
in classification - and then next, we'll switch to the midpoint filtering
just for the sake of getting numbers

that will then need to be repeated with nuclear filtering to see the difference
made, after which I iteratively add more and more filters for a sense of how
much is retained/lost 

on top of this, I need to get the numbers on how much was lost at each
point from a given filter - fortunately this is what all the logs 
were implemented for! going to go through these in a separate analysis
folder called `filters` 

in the meantime, going to

- remove mapq and base qual filter in `filter` (forgot I had the latter here!) 
- switch to midpoint filtering
- turn off base quality filter in `classification` (via `summarise_cross`)
- turn off base quality filter in `parental_phase_changes.smk`? 

going to do this via `phase_change_detection.smk` - going to leave
`phase_change_nuclear.smk` as is for now

finally, need to update `summarise_cross` so that it also outputs
base qualities of phase change containing variants

```bash
rm -v data/phase_changes/event_summaries
rm -v data/phase_changes/sam
rm -v data/phase_changes/bam
rm -v data/phase_changes/fp_bed

```

changes to readcomb:
- add base qual filter to command line args in readcomb-fp
- add base qual as attribute to `Pair` objects 
    - mean base qual and min base qual for phase-change-involved variants

will need to see the distribution of base quals after the fact as well - but let's
get re-filtering first, with all min mapq and base qual filters set to 0

going to just get readcomb-filter going again while I rejig readcomb-fp and classification:

```bash
# rule all just set to the 'init.sam' files, no fp filtering here
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 24
```

and now to get these readcomb changes going - I've updated the detection
method to include base quals, switched the default base qual filter to 0, and
added `min_base_qual`, `mean_base_qual`, and `phase_change_variants` as new
attributes to `Pair` objects

recompiling and then off we go - the workflow should continue running as expected

in the meantime, I'm going to work on `summarise_cross` - let me see if I can
implement multiprocessing here to speed it up somewhat

alright, blind run after a test implementation:

```bash
mkdir p_test
head -n 10000 data/phase_changes/sam/2344x1952.init.sam > p_test/input.sam

python analysis/phase_changes/summarise_cross.py \
--bam p_test/input.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--mask_size 75 --processes 10 --mapq 0 --base_qual 20 \
--remove_uninformative --out p_test/csv_test
```

after an hour of trial and error and various bug fixes... this works!!! 

going to update it tomorrow to handle and concatenate the temp files but I'm so excited!

## 10/8/2022

today - 

- updating `summarise_cross` to concatenate and handle files correctly
- if possible - implementing a counter queue for `summarise_cross`
- doing a speed test on `readcomb-fp` and then parallelizing if necessary

the first step should be pretty straightforward - going to start with that and then test

here goes:

```bash
time python analysis/phase_changes/summarise_cross.py \
--bam p_test/input.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--mask_size 75 --processes 4 --mapq 0 --base_qual 20 \
--remove_uninformative --out p_test/csv_test
```

looks great! now for a quick implementation of a counter object and then we're good to go

next up - a speed test on readcomb-fp (although these fp files will need to be
regenerated downstream I think)

```bash
time readcomb-fp --fname data/phase_changes/sam/2344x1952.init.sam \
--false_plus data/phase_changes/parental/2344x1952.plus.filtered.sam \
--false_minus data/phase_changes/parental/2344x1952.minus.filtered.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--method midpoint --read_length_filter 0 --log fp.log \
--false_bed_out 2344x1952.test.bed.out --out 2344x1952.test.filter.sam
# took 36 min - 6 min for bed creation, 30 min for actual filtering
```

seems to be a bug here - need to update the tabix command to include `-b2 -e3`

## 11/8/2022

going to recompile readcomb with that bug fix and then get the fp stuff running 

`readcomb-filter` completed in 18 hours (instead of the 15 days it would have taken otherwise!)
so that's pretty good - just need fp done by tonight, after which I'll run the updated
`summarise_cross` and get some descriptive stats about quality metrics and the distributions
of metrics like mapq

though it occurs to me I still need to re-run readcomb-filter on the parental files! let's
get that going -

```bash
time snakemake -pr -s analysis/phase_changes/parental_phase_changes.smk --cores 30
# going very hard on the server here but I'm in a rush and it seems clear!  
# min mapq and base qual set to 0 - although I should get the mapq distribution of parental reads
```

getting the mapq distribution (among other things) for these parental sequences - 
gonna get a Jupyter notebook going for this - and a new folder called `filters`
to get these sorts of statistics in

## 12/8/2022

lol that took 27 hours

now for the actual fp step:

```bash
# rule all now set to actually do fp filtering
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 8
# took 18 hours - and that includes the sorting/bam conversion/indexing after
```

will also create bams for IGV viewing if needed:

```bash
# in data/phase_changes/
for fname in parental/*sam; do
    f=$(basename ${fname} .sam)
    samtools sort -O bam -o parental/${f}.bam ${fname};
    samtools index parental/${f}.bam;
done
```

## 13/8/2022

today - 

- run `summarise_cross` on all samples
- get metrics of reads lost along the various filtering stages (e.g. by looking at filter/fp logs)
- generate single table with summaries of this - will probably do this in an Rmd
- requeue nuclear filtering in `phase_changes/nuclear` - this is to contrast nuclear/midpoint filtering

running `summarise_cross` - now supports multiprocessing, so that needs to be added as
a parameter in the workflow

```bash
# rule all now updated to generate tsvs
# no longer has --remove_uninformative either for easier overall summary stat calculations
# had to upload readcomb to handle reads that are too short vs masking size correctly

time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

all done, and the parallelization meant that it just took 11 minutes to go through
all 32 crosses! not bad at all

redoing the nuclear filtering as well to contrast the numbers:

```bash
# modified the workflow to use the init files in data/phase_changes/sam
# since these should be identical either way - it's just the fp method that changes

time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 12
```

next up worth investigating - what is up with read pairs that `filter` considers
phase changes but which are later deemed `no_phase_change`? what causes this? 

```python
import readcomb.classification as rc
from tqdm import tqdm

bam_fname = 'data/phase_changes/sam/2344x1952.init.sam'
vcf_fname = 'data/genotyping/vcf_filtered/2344x1952.vcf.gz'

reader = rc.pairs_creation(bam_fname, vcf_fname)

weird = []
for pair in tqdm(reader):
    pair.classify(masking=0, quality=0)
    if pair.call == 'no_phase_change':
        weird.append(pair)
    if len(weird) == 10:
        break

x = weird[0]

>>> print(x)
Record name: A00516:220:H2NGGDRXY:1:2101:1217:32565
Read1: chromosome_16:1090803-1090877
Read2: chromosome_16:1090826-1091076
VCF: data/genotyping/vcf_filtered/2344x1952.vcf.gz
Unmatched Variant(s): False
Condensed: [['CC1952', 1090803, 1091076]]
Call: no_phase_change
Condensed Masked: [['CC1952', 1090803, 1091076]]
Masked Call: no_phase_change
Midpoint: ('chromosome_16:1090939', -1)
Variants Per Haplotype: 9.0
Variant Skew: -1
Gene Conversion Length: N/A
Min Variants In Haplotype: -1
Mismatch/Variant Ratio: 1.7777777777777777

>>> x.detection
[('2', 1090814, 'A', 37), ('2', 1090825, 'G', 37), ('2', 1090835, 'C', 37), 
('2', 1090938, 'T', 11), ('2', 1090940, 'A', 37), ('2', 1091003, 'T', 37), 
('2', 1091062, 'C', 37), ('2', 1091064, 'T', 25), ('2', 1091065, 'G', 37)]

# and now to check against the readcomb filter algorithm
from readcomb.filter import cigar, qualities_cigar, phase_detection
import argparse

parser = argparse.ArgumentParser
parser.quality = 0

>>> phase_detection(x.variants_1, cigar(x.rec_1), x.rec_1, parser)
['2', '2', '2', '2']

>>> phase_detection(x.variants_2, cigar(x.rec_2), x.rec_2, parser)
['2', '1', '2', '2', '2', '2', '2', '2'] # weird! 

# looks like was an overlap issue
>>> rc.downstream_phase_detection(x.variants_1, cigar(x.rec_1), x.rec_1, 0)
[('2', 1090814, 'A', 37), ('2', 1090825, 'G', 37), ('2', 1090835, 'C', 37), ('2', 1090840, 'G', 37)]

>>> rc.downstream_phase_detection(x.variants_2, cigar(x.rec_2), x.rec_2, 0)
[('2', 1090835, 'C', 37), ('1', 1090840, 'C', 11), ('2', 1090938, 'T', 11), 
('2', 1090940, 'A', 37), ('2', 1091003, 'T', 37), ('2', 1091062, 'C', 37), 
('2', 1091064, 'T', 25), ('2', 1091065, 'G', 37)]
# see site 1090840

# is this consistently the case? 
x = weird[1]
>>> rc.downstream_phase_detection(x.variants_1, cigar(x.rec_1), x.rec_1, 0)
[('2', 4393536, 'T', 25), ('2', 4393542, 'T', 37), ('2', 4393575, 'G', 37), 
('2', 4393576, 'C', 37), ('2', 4393582, 'G', 37), ('2', 4393620, 'A', 37), 
('2', 4393636, 'G', 37), ('2', 4393640, 'G', 37), ('2', 4393657, 'C', 37), 
('2', 4393678, 'G', 37), ('2', 4393688, 'A', 37), ('2', 4393717, 'G', 25), ('2', 4393718, 'A', 37)]
>>> rc.downstream_phase_detection(x.variants_2, cigar(x.rec_2), x.rec_2, 0)
# see last record on first line
[('2', 4393536, 'T', 37), ('2', 4393542, 'T', 37), ('2', 4393575, 'G', 37), ('1', 4393576, 'G', 11), 
('2', 4393582, 'G', 37), ('2', 4393620, 'A', 37), ('2', 4393636, 'G', 37), ('2', 4393640, 'G', 37), 
('2', 4393657, 'C', 37), ('2', 4393678, 'G', 37), ('2', 4393688, 'A', 37), ('2', 4393717, 'G', 37), 
('2', 4393718, 'A', 25)]

# let's iterate through the file and see if there's a `no_phase_change` outside of an overlap read
reader = rc.pairs_creation(bam_fname, vcf_fname)
more_weird = []
for pair in tqdm(reader):
    pair.classify(masking=0, quality=0)
    if pair.call == 'no_phase_change' and not pair.overlap:
        more_weird.append(pair)
        print('got one')
    if len(more_weird) == 10:
        break # took 8250 iterations to find 10

# seems all of these end in a single N variant - which is bugging out the creation of of Pair.condensed
# let's do a quick check to see whether Ns and overlaps explain it all:

most_weird = []
reader = rc.pairs_creation(bam_fname, vcf_fname)
for pair in tqdm(reader):
    pair.classify(masking=0, quality=0)
    if pair.call == 'no_phase_change' and not pair.overlap and not pair.no_match:
        most_weird.append(pair)
        print('got one')
    if len(most_weird) == 10:
        break # took 107000 iterations over 14 min to find 10
```

## 14/8/2022

figured it out - `Pair._create_condensed` needs to have a separate
copy of `detection` that ignores all Ns (explains `more_weird`) and
I need to add a stipulation that the 'false haps' where start = end in 
detection are only removed if the reads overlap (since otherwise 
single base pair haps caused by the read starting with two variants
are removed)

redoing readcomb-fp, because of course - 

```bash
# going to move the previous files elsewhere for now
mkdir -p data/phase_changes/midpoint_temp/
mv -v data/phase_changes/sam/*filtered.sam data/phase_changes/midpoint_temp/
mv -v data/phase_changes/fp_bed/* data/phase_changes/midpoint_temp/
mv -v data/phase_changes/event_summaries/* data/phase_changes/midpoint_temp/

rm -v data/phase_changes/bam/*

time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 8
```

and same goes for the nuclear filtering stuff:

```bash
# in data/phase_changes/nuclear
mkdir nuclear_temp
mv -v sam/*filtered.sam nuclear_temp/
mv -v fp_bed/* nuclear_temp/
mv -v event_summaries/* nuclear_temp/

rm -v bam/*

# back to root
time snakemake -pr -s analysis/phase_changes/phase_change_nuclear.smk --cores 8
```

## 15/8/2022

both done in 2 hours! not so bad - not sure why I didn't run these in parallel
on snakemake sooner

exporting the logs to compare midpoint vs nuclear filtering, as well as the output
of `summarise_cross`, in an Rmd

## 16/8/2022

today - going to run `summarise_cross` on the parents as well to make
sense of why some have more false phase changes than others:

```bash
# new rule that generates these in data/phase_changes/parental_event_summaries
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

## 18/8/2022

today - need to implement one extra quality metric that specifies whether
a read would have been removed by the nuclear filter or not (e.g. if it 
overlaps a false phase change read at all, even if the phase change itself doesn't
actually match - since we're sticking to the midpoint option) 

after that, going to go back to labelling events to then train a model on

quick test:

```bash
python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1691.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1691.nuclear.bed.gz \
--mask_size 75 --processes 12 --mapq 0 --base_qual 0 \
--out sc_test.tsv
```

looks good - let's clear out `data/phase_changes/event_summaries`, update the workflow,
and regen these files

```bash
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

## 24/8/2022

coming back to this just to add a `min_phase_change_var_qual` attribute to `summarise_cross`

testing:

```bash
python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1691.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1691.nuclear.bed.gz \
--mask_size 75 --processes 12 --mapq 0 --base_qual 0 \
--out sc_test.tsv
```

looks good - let's do this again -

```bash
# clearing old ones and doing this again
rm -v data/phase_changes/event_summaries/
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

and then going to take this offline so I can use my labels with it

## 29/8/2022

today - checking to see whether a suspected bug in `summarise_cross` is as problematic
as I think it is - it seems the last pair is being ignored, which means if there were
n processes that n pairs have been excluded

let's just run an example but with 12 processes and then 16 - 

```bash
python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1691.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1691.nuclear.bed.gz \
--mask_size 75 --processes 12 --mapq 0 --base_qual 0 \
--out sc_test_12.tsv

python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1691.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1691.nuclear.bed.gz \
--mask_size 75 --processes 16 --mapq 0 --base_qual 0 \
--out sc_test_16.tsv

# if this bug is legit these should have different line counts
wc -l sc_test_12.tsv # 34129
wc -l sc_test_16.tsv # 34129
# looks alright! 
```

next up - trying to figure out why some crossovers are being ignored by `summarise_cross` -
I had a previously labelled set of crossovers from 2343x1952 but they seem to be missing

```bash
python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1952.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1952.nuclear.bed.gz \
--mask_size 75 --processes 12 --mapq 0 --base_qual 0 \
--out sc_test.tsv
```

what on earth?

```
$ wc -l sc_test.tsv
36994 sc_test.tsv

hasans11@hpcnode1 /research/projects/chlamydomonas/genomewide_recombination/rcmb (work)
$ wc -l data/phase_changes/event_summaries/2343x1952.75.tsv
36994 data/phase_changes/event_summaries/2343x1952.75.tsv

hasans11@hpcnode1 /research/projects/chlamydomonas/genomewide_recombination/rcmb (work)
$ grep -c 'A00516:220:H2NGGDRXY:2:2120:4029:9377' sc_test.tsv
1

hasans11@hpcnode1 /research/projects/chlamydomonas/genomewide_recombination/rcmb (work)
$ grep -c 'A00516:220:H2NGGDRXY:2:2120:4029:9377' data/phase_changes/event_summaries/2343x1952.75.tsv
0
```

same number of events? but only one has the read of interest? 

I have no idea what's going on - but let's delete this file and get snakemake to regen it
before testing again

```bash
hasans11@hpcnode1 /research/projects/chlamydomonas/genomewide_recombination/rcmb (work)
$ grep -c 'A00516:220:H2NGGDRXY:2:2102:13928:3317' data/phase_changes/event_summaries/2343x1952.75.tsv
1

# looks good now??
```

let's just rerun all of these to be on the safe side

```bash
rm -v data/phase_changes/event_summaries/*tsv
snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

## 3/9/2022

today - trying to get tandem repeats finder working

looks like I need parental fastas for this, which were otherwise generated in
`analysis/ldhelmet` - let's give this a go for a file that contains one
chrom across all parents:

```bash
time ./bin/trf data/ldhelmet/fasta_chrom/NA1.chromosome_01.fa \
2 7 7 80 10 500 -d
```

alright! so this keeps doing this thing where it saves repeats in intermediate files,
deletes them all at the end, and then proudly proclaims there are no repeats

going to switch over to a tool called ULTRA - will need to mount a docker image for this:

```bash
mkdir analysis/tandem_repeats/
# in the above folder - using the globally installed singularity exec
singularity build ultra docker://traviswheelerlab/ultra
singularity run ultra_latest.sif # looks good
```

testing on a fasta -

```bash
# need absolute paths for this to work
time singularity exec --bind /research/projects/chlamydomonas/genomewide_recombination/rcmb/data/ldhelmet/fasta_chrom/:\
/research/projects/chlamydomonas/genomewide_recombination/rcmb/analysis/tandem_repeats/mnt \
analysis/tandem_repeats/ultra.sif \
ultra -f test.json -at 0.36 \ # GC content
-s 9 data/ldhelmet/fasta_chrom/NA2.chromosome_01.fa # -s is score threshold
```

## 4/9/2022

looks good - checking to see whether I can parse this in python - 

```python
import json
with open('test.json', 'r') as f:
    reader = json.load(f)
    x = reader['Repeats'][0]
print(x) # prints first line as dict
```

this can be pretty easily converted to a tsv then -

```python
import json
import csv
from tqdm import tqdm

with open('test.json', 'r') as f_in:
    with open('test_converted.tsv', 'r') as f_out:
        reader = json.load(f_in)
        fieldnames = list(reader['Repeats'][0].keys())
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        
        for record in tqdm(reader['Repeats']):
            writer.writerow(record)
```

this looks good - but I need to convert these into tabix-able files for `summarise_cross` to
get metrics for each pair as needed

going to continue this in `analysis/tandem_repeats/log.md`

--
            
alright - now that these have been generated, I've updated `summarise_cross` to return
various microsat metrics

test run -

```bash
time python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1952.vcf.gz \
--false_bed data/phase_changes/nuclear/fp_bed/2343x1952.nuclear.bed.gz \
--microsats data/tandem_repeats/samples/CC2343.sorted.tsv.gz \
data/tandem_repeats/samples/CC1952.sorted.tsv.gz \
--mask_size 75 --processes 12 --mapq 0 --base_qual 0 \
--out sc_test.tsv
```

took a bit of debugging but looks good - just need to update the workflow
to correctly pull the microsat files

```bash
# this time with microsats
time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk --cores 12
```

it does occur to me that the parental event summaries can't be run again without
microsats being implemented as an optional arg - but I'll worry about that if/when I have to

## 5/9/2022

quickly checking to see what's up with the index files - let's just do a test run:

```bash
# going to cancel out soon as the bed file is made
time readcomb-fp --fname data/phase_changes/sam/2344x1952.init.sam \
--false_plus data/phase_changes/parental/2344x1952.plus.filtered.sam \
--false_minus data/phase_changes/parental/2344x1952.minus.filtered.sam \
--vcf data/genotyping/vcf_filtered/2344x1952.vcf.gz \
--method nuclear --read_length_filter 0 --log fp.log \
--false_bed_out fp.temp.bed --out fp.temp.sam
```

## 12/9/2022

hello from Northern ON!

index files were fine - I did a diff that day and there was nothing
different about them

today - reworked min end proximity to handle overlapping sites better.
this is because this quality metric was lower than it should
have been for phase change variants found in an overlapping region

going to do a test run on 2343x1952 since that had my test case
read

```bash
time python analysis/phase_changes/summarise_cross.py \
--bam data/phase_changes/sam/2343x1952.filtered.sam \
--vcf data/genotyping/vcf_filtered/2343x1952.vcf.gz \
--mask_size 75 --processes 16 --mapq 0 --base_qual 0 \
--microsats data/tandem_repeats/samples/CC2343.sorted.tsv.gz \
data/tandem_repeats/samples/CC1952.sorted.tsv.gz \
--out sc_test.tsv
```

tested for differences by loading both this and the old file into R
and seeing whether the min end prox column vectors differed on average
(eg mean difference of both vectors) - need to sort files first though
of course, since they were generated w/ multiprocessing

```R
library(tidyverse)
d_old = read_tsv('data/phase_changes/event_summaries/2343x1952.75.tsv', 
    col_types = cols()) %>%
    arrange(midpoint, read_name)
d_new = read_tsv('sc_test.tsv') %>%
    arrange(midpoint, read_name)
```

looks good - this only seems to affect 6197 of the 73989 rows

going to rerun this across all crosses and then redo the model tomorrow

```bash
rm -v data/phase_changes/event_summaries/*

time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk \
--cores 16
```

## 13/9/2022

alright - exporting these offline and rerunning the model

next up - redoing the microsat files with much lower thresholds - the model
should be able to sort through 'low quality' microsats and whether they
correlate with false calls or not - I just need to give it more data to work
with! 

going to quickly do this by rerunning the tandem repeats calc workflow
but with the score threshold set to 0

```bash
rm -v data/tandem_repeats/json/*
rm -v data/tandem_repeats/samples/*
rm -v data/tandem_repeats/tsv/*

time snakemake -pr -s analysis/tandem_repeats/tandem_repeats_calc.smk \
--cores 8
```

and now to regen the event summaries - again:

```bash
rm -v data/phase_changes/event_summaries/*

time snakemake -pr -s analysis/phase_changes/phase_change_detection.smk \
--cores 16
```

## 20/9/2022

today - trying to fit a random forest regression tree to the training data -
I would do this locally but I can't get the `randomForest` package working on my local
machine, so I'm going to try this in a Jupyter Notebook on here

alright nevermind - Jupyter is completely shitting the bed cause of some X11 issue
I don't have the brainpower or desire to make sense of, so we're doing this in a kernel

```R
library(tidyverse)
library(tree) # had to install.packages
library(randomForest)
set.seed(42)
d = read_tsv('data/phase_changes/cos_all_annotated.tsv', col_types = cols())
model_vars = d %>%
    select(-comments, -cross, -chromosome, -midpoint, -call,
        -start, -end, -mask_size, -masked_call, -detection, -gc_length) %>%
    mutate(
        check = factor(check),
        outer_bound = ifelse(outer_bound > 0.5, 1 - outer_bound, outer_bound),
        rel_midpoint = ifelse(rel_midpoint > 0.5, 1 - rel_midpoint, rel_midpoint),
        indel_proximity = ifelse(indel_proximity == -1, pmax(read1_length, read2_length), indel_proximity),
        proximate_indel_length = ifelse(proximate_indel_length == -1, 0, proximate_indel_length),
        false_overlap = factor(false_overlap),
        indel_close = ifelse(indel_proximity <= 5, 1, 0),
        parent1_microsat = factor(parent1_microsat, levels = c(0, 1)),
        parent2_microsat = factor(parent2_microsat, levels = c(0, 1)),
        parent1_microsat_score = ifelse(is.na(parent1_microsat_score), 0, parent1_microsat_score),
        parent2_microsat_score = ifelse(is.na(parent2_microsat_score), 0, parent2_microsat_score),
        parent1_microsat_log2_pval = ifelse(is.na(parent1_microsat_log2_pval), 0, parent1_microsat_log2_pval),
        parent2_microsat_log2_pval = ifelse(is.na(parent2_microsat_log2_pval), 0, parent2_microsat_log2_pval),
        parent1_microsat_proximity = ifelse(
          is.na(parent1_microsat_proximity), pmax(read1_length, read2_length), parent1_microsat_proximity),
        parent2_microsat_proximity = ifelse(
          is.na(parent2_microsat_proximity), pmax(read1_length, read2_length), parent2_microsat_proximity),
        min_base_qual_discrete = factor(
          min_base_qual, levels = c('11', '25', '37'),
          labels = c('11', '25', '37'), ordered = TRUE)
     )

test_subset = model_vars %>% slice_sample(n = 105) 
train_subset = model_vars %>% anti_join(test_subset, by = 'read_name')

test_subset = test_subset %>% select(-read_name)
train_subset = train_subset %>% select(-read_name)

rf_model = randomForest(check ~ ., data = train_subset, importance = TRUE, proximity = TRUE)
print(rf_model)

Call:
 randomForest(formula = check ~ ., data = train_subset, importance = TRUE,      proximity = TRUE)
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 6

        OOB estimate of  error rate: 6.14%
Confusion matrix:
    0   1 class.error
0 783  10  0.01261034
1  48 104  0.31578947

x = predict(rf_model, newdata = test_subset)
y = test_subset$check
table(x, y)
   y
x    0  1
  0 89  7
  1  2  7
# not great! 

# training error rate?
rf_training = randomForest(check ~ ., data = model_vars, importance = TRUE, proximity = TRUE)
Call:
 randomForest(formula = check ~ ., data = model_vars, importance = TRUE, proximity = TRUE)
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 6

        OOB estimate of  error rate: 6.76%
Confusion matrix:
    0   1 class.error
0 866  18  0.02036199
1  53 113  0.31927711

```

## 25/9/2022

today - need to get effective sequence per window across the 
genome for each of the crosses, since this will give me the denominator
when calculating windowed CO rates 

one way to do this - calculate midpoints independently from
readcomb as just `(end - start) / 2` while otherwise getting
end and start the same way as I did when generating the `read_lengths.tsv`
file earlier (9/7/2022)

but one problem is that these BAM files aren't sorted... 
let's try reimplementing that 17 array set idea that I had

better yet, could create and repeatedly update a single giant 19 x n array,
where n is the maximum number of possible windows - this avoids having to deal
with a Python dict

could do two arrays so I could actually keep the number of read counts as well
come to think of it

one concern is gene conversions - but I could always subtract these from the
values after the fact

going to get started on `windowed_effective_sequence.py`

giving this a go:

```bash
time python analysis/phase_changes/windowed_effective_sequence.py \
--fname data/alignments/bam_prepped/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--window_size 10000 \
--out eff_test.tsv
```

## 26/9/2022

scaling this up to all crosses - although I think I should
start a new analysis folder called `rate` for this, and I'm going
to move this script there in anticipation of more landscape data

## 27/9/2022

need to do LOOCV cross validation, but it's probably way too intensive to
run on my local machine so I'm going to do it here

```R
library(plyr)
library(tidyverse)
library(glmnet)

all_annotated = read_tsv('data/phase_changes/cos_all_annotated.tsv', col_types = cols()) %>%
     mutate(
    outer_bound = ifelse(outer_bound > 0.5, 1 - outer_bound, outer_bound),
    rel_midpoint = ifelse(rel_midpoint > 0.5, 1 - rel_midpoint, rel_midpoint),
    indel_proximity = ifelse(indel_proximity == -1, pmax(read1_length, read2_length), indel_proximity),
    proximate_indel_length = ifelse(proximate_indel_length == -1, 0, proximate_indel_length),
    false_overlap = factor(false_overlap),
    indel_close = ifelse(indel_proximity <= 5, 1, 0),
    parent1_microsat = factor(parent1_microsat, levels = c(0, 1)),
    parent2_microsat = factor(parent2_microsat, levels = c(0, 1)),
    parent1_microsat_score = ifelse(is.na(parent1_microsat_score), 0, parent1_microsat_score),
    parent2_microsat_score = ifelse(is.na(parent2_microsat_score), 0, parent2_microsat_score),
    parent1_microsat_log2_pval = ifelse(is.na(parent1_microsat_log2_pval), 0, parent1_microsat_log2_pval),
    parent2_microsat_log2_pval = ifelse(is.na(parent2_microsat_log2_pval), 0, parent2_microsat_log2_pval),
    parent1_microsat_proximity = ifelse(
      is.na(parent1_microsat_proximity), pmax(read1_length, read2_length), parent1_microsat_proximity),
    parent2_microsat_proximity = ifelse(
      is.na(parent2_microsat_proximity), pmax(read1_length, read2_length), parent2_microsat_proximity),
    min_base_qual_discrete = factor(
      min_base_qual, levels = c('11', '25', '37'),
      labels = c('11', '25', '37'), ordered = TRUE),
    check = factor(check)) 

# initial lasso fit
set.seed(42)
model_vars = all_annotated %>% 
  select(-comments, -cross, -read_name, -chromosome, -midpoint, -call,
         -start, -end, -mask_size, -masked_call, -detection, -gc_length)
  # mutate(check = as.numeric(check))
pred_matrix = model.matrix(check ~ .*., model_vars)
outcome = model_vars$check

# obtain lambda value
cv_lasso = cv.glmnet(
  pred_matrix, outcome, type.measure = 'deviance', alpha = 1, family = 'binomial')

lasso_model = glmnet(
  pred_matrix, outcome, family = 'binomial', 
  alpha = 1, # lasso regression
  lambda = cv_lasso$lambda)

# get training probs - supply pred matrix
lasso_training_probs = predict(
  lasso_model, newx = pred_matrix, s = cv_lasso$lambda.min, type = 'response')

# create vector of predictions
lasso_training_preds = ifelse(lasso_training_probs > 0.4, 1, 0)

# confusion matrix
table(lasso_training_preds, all_annotated$check)

mean(lasso_training_preds == all_annotated$check)

# LOOCV
cv_preds = numeric(length = nrow(all_annotated))
cv_probs = numeric(length = nrow(all_annotated))

pbar = create_progress_bar('text')
pbar$init(nrow(all_annotated))

for (i in 1:nrow(all_annotated)) {
    train_split = all_annotated[-i,]
    test_split = all_annotated[i,]
    model_vars = train_split %>%
        select(-read_name, -chromosome, -midpoint, -call, -start, -end,
               -comments, -cross, -mask_size, -masked_call, -detection, -gc_length)
    pred_matrix = model.matrix(check ~ .*., model_vars)

    cv_lasso_current = cv.glmnet(
        pred_matrix, model_vars$check, type.measure = 'deviance', alpha = 1, family = 'binomial')
    lasso_model_current = glmnet(
        pred_matrix, model_vars$check, family = 'binomial',
        alpha = 1, lambda = cv_lasso_current$lambda)
    
    test_model_vars = test_split %>%
        select(-read_name, -chromosome, -midpoint, -call, -start, -end,
               -comments, -cross, -mask_size, -masked_call, -detection, -gc_length)
    test_pred_matrix = model.matrix(check ~ .*., test_model_vars)
    cv_probs[i] = predict(
        lasso_model_current, newx = test_pred_matrix, 
        s = cv_lasso_current$lambda.min, type = 'response')
    cv_preds[i] = ifelse(cv_probs[i] > 0.4, 1, 0)
    pbar$step()
}

table(cv_preds, all_annotated$check, deparse.level = 2)

        # all_annotated$check
# cv_preds   0   1
       # 0 849  32
       # 1  35 134

all_annotated$loocv_pred = cv_pred
all_annotated$loocv_prob = cv_probs

write_tsv(all_annotated, 'data/phase_changes/cos_all_annotated_loocv.tsv')

```
