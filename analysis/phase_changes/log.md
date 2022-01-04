
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
