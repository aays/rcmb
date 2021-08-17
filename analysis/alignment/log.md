
## 29/6/2021

time to download those sweet, sweet fastq files - doing this in `scratch` for now
given space issues

see Notion for details - left running overnight

## 30/6/2021

files are now in `/scratch/projects/chlamydomonas/recombination/fastq/`

creating symlinks:

```bash
mkdir -p data/alignments/fastq
ln -sv /scratch/projects/chlamydomonas/recombination/fastq/* .
mv -v ../readSets.md5 # brought over from local machine, downloaded from Nanuq
```

and now for the checksum comparison:

```bash
time md5sum -c readSets.md5
```

looks good:

```
NS.A00516_0220.IDT_i7_47---IDT_i5_47.2932x2342_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_47---IDT_i5_47.2932x2342_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_9---IDT_i5_9.2932x1691_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_9---IDT_i5_9.2932x1691_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_21---IDT_i5_21.2344x3059_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_21---IDT_i5_21.2344x3059_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_157---IDT_i5_157.2932x2935_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_157---IDT_i5_157.2932x2935_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_121---IDT_i5_121.3086x1691_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_121---IDT_i5_121.3086x1691_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_45---IDT_i5_45.2343x3059_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_45---IDT_i5_45.2343x3059_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_145---IDT_i5_145.GB119x2342_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_145---IDT_i5_145.GB119x2342_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_10---IDT_i5_10.3071x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_10---IDT_i5_10.3071x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_33---IDT_i5_33.2343x2935_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_33---IDT_i5_33.2343x2935_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_93---IDT_i5_93.3071x2342_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_93---IDT_i5_93.3071x2342_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_12---IDT_i5_12.2344x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_12---IDT_i5_12.2344x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_35---IDT_i5_35.2932x2931_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_35---IDT_i5_35.2932x2931_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_94---IDT_i5_94.GB119x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_94---IDT_i5_94.GB119x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_24---IDT_i5_24.2343x1691_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_24---IDT_i5_24.2343x1691_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_133---IDT_i5_133.3086x3059_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_133---IDT_i5_133.3086x3059_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_71---IDT_i5_71.2344x3062_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_71---IDT_i5_71.2344x3062_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_95---IDT_i5_95.2344x2342_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_95---IDT_i5_95.2344x2342_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_57---IDT_i5_57.2343x2931_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_57---IDT_i5_57.2343x2931_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_23---IDT_i5_23.2932x3062_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_23---IDT_i5_23.2932x3062_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_46---IDT_i5_46.3086x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_46---IDT_i5_46.3086x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_59---IDT_i5_59.2932x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_59---IDT_i5_59.2932x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_81---IDT_i5_81.3071x3059_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_81---IDT_i5_81.3071x3059_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_22---IDT_i5_22.3086x2935_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_22---IDT_i5_22.3086x2935_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_70---IDT_i5_70.GB119x2935_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_70---IDT_i5_70.GB119x2935_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_58---IDT_i5_58.GB119x1691_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_58---IDT_i5_58.GB119x1691_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_36---IDT_i5_36.2343x2342_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_36---IDT_i5_36.2343x2342_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_82---IDT_i5_82.GB119x3062_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_82---IDT_i5_82.GB119x3062_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_34---IDT_i5_34.3086x2931_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_34---IDT_i5_34.3086x2931_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_69---IDT_i5_69.2343x1952_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_69---IDT_i5_69.2343x1952_R2.fastq.gz: OK
NS.A00516_0220.IDT_i7_83---IDT_i5_83.2344x2931_R1.fastq.gz: OK
NS.A00516_0220.IDT_i7_83---IDT_i5_83.2344x2931_R2.fastq.gz: OK
```

## 14/7/2021

finally getting back on this - time to rename these symlinks to make life a bit easier

```python
import os
import re
from glob import glob

# in data/alignments/fastq
fq_list = glob.glob('*R?.fastq.gz')

pattern = '[A-Z0-9]{3,}x[A-Z0-9]{3,}_R[12].fastq.gz'
for fname in fq_list:
    new_fname = re.search(pattern, fname).group()
    os.rename(fname, new_fname)
```

looks good - now to get snakemake going

## 16/7/2021

first - getting references and the dir in order

```bash
bash analysis/alignment/dir_setup.sh

# in data/references
ln -sv /research/references/chlamydomonas/6.0_chlamy/* .
```

next up - fastqc:

```bash
conda install -c bioconda fastqc
time snakemake -p --snakefile analysis/alignment/read_qc_trim.smk
```

## 17/7/2021

ran for 7 hours but can't argue with results! 

now to unzip the raw fastqc zip files - gotta be careful to
avoid fastqc starting up again

```bash
time snakemake -p \
--snakefile analysis/alignment/read_qc_trim.smk \
unzip_fastqc
```

looking for failed tests:

```bash
cd data/alignments/fastq/fastqc/
grep -P 'FAIL\t' */summary.txt
```

looks like all of thse are adapter content or GC content failures - let's
look for per base sequence issues explicitly:

```bash
grep -P 'FAIL\tPer base sequence quality' */summary.txt # nothing! 
```

next up - will have to per Trimmomatic, get the adapters, and trim reads

## 19/7/2021

adapters are in `/research/tmp_apps/Trimmomatic/adapters/` - specifically
`NEBNext_dual.fasta`

```python
>>> from Bio import SeqIO
>>> d = [s for s in SeqIO.parse('/research/tmp_apps/Trimmomatic-0.36/adapters/NEBNext_dual.fasta', 'fasta
... ')]
>>> [s.id for s in d]
['NEBNext_dualRead_1', 'NEBNext_dualRead_2', 'NEBNext_dualRead_1_RC', 'NEBNext_dualRead_2_RC']

>>> x = ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT']
>>> x[0] in [str(s.seq) for s in d]
True

>>> x[1] in [str(s.seq) for s in d]
True
```

also includes the reverse complements, which is nice! symlinked into `bin`
and it's off to the races w/ snakemake

```bash
snakemake -p --snakefile analysis/alignment/read_qc_trim.smk \
--cores 16 trim_reads # going to take a while, so using 16 cores
```

## 21/7/2021

trimming finally done - took a day and 6 hours, and that's with
16 threads! 

getting fastqc going again:

```bash
time snakemake -p --snakefile analysis/alignment/read_qc_trim.smk fastqc_trim
```

in the meantime, used tmux to copy stderr from trimmomatic into
`data/alignments/fastq_trim/trim.log` - forgot to set up the redirect! 

## 23/7/2021

fastqc'ing took another 13 hours but done now - time to check again:

```bash
time snakemake -p \
--snakefile analysis/alignment/read_qc_trim.smk \
unzip_fastqc_trim
```

looking for failed tests this time around:

```bash
cd data/alignments/fastq_trim/fastqc/
grep -P 'FAIL\t' */summary.txt
```

only expected GC content 'failures' and a few for per base seq content

making sure:

```bash
grep -P 'FAIL\tPer base' */summary.txt # no issues with sequence quality! 
```

looks good - time to align! starting a new snakefile called `alignment.smk`

going to make a quick list of filenames to expedite snakemake rule writing:

```bash
touch data/alignments/fastq_trim/trim_names.txt
touch data/alignments/samples.txt

for fname in data/alignments/fastq_trim/*x*_trim_?.fq.gz; do
    basename ${fname .fq.gz} >> data/alignments/fastq_trim/trim_names.txt;
done

for fname in data/alignments/fastq_trim/*x*_trim_?.fq.gz; do
    basename ${fname%_trim_*.fq.gz} >> data/alignments/samples.txt;
done

sort data/alignments/samples.txt | uniq > data/alignments/samples_final.txt
rm data/alignments/samples.txt
mv -v data/alignments/samples_final.txt data/alignments/samples.txt

# also this
mkdir -p data/alignments/bam_temp
```
    
and now to get this show on the road:

```bash
time snakemake -p -s analysis/alignment/alignment.smk --cores 16
```

## 26/7/2021

took 22 hours but we're in the clear! 

now to sort the bams:

```bash
# rule bam_sort
time snakemake -pr -s analysis/alignment/alignment.smk --cores 4
```

fix mate info:

```bash
# rule bam_fix_mate
time snakemake -pr -s analysis/alignment/alignment.smk # took 15 hrs
```

## 27/7/2021

add read groups - each sample had two libraries, so we don't need to
be too worried about lane effects

```bash
time snakemake -pr -s analysis/alignment/alignment.smk # done in 7 hrs
```

## 28/7/2021

time for MarkDuplicates -

```bash
time snakemake -pr -s analysis/alignment/alignment.smk
``` 

## 9/8/2021

so I assumed this had completed and moved on to salt project stuff for a bit,
but it never actually did since it ran out of disk space! `temp()` in the snake
file wasn't clearing the intermediate files in `bam_temp` (likely cause I was
writing the rules one at a time) so I had to manually remove them just now

here goes once more:

```bash
time snakemake -pr -s analysis/alignment/alignment.smk # took 16 hours
```

this should finish just fine - but I still need to realign the parental
sequences to v6 for variant calling as well! can reuse much
of `alignment.smk` but going to make a copy anyways so
any other changes can be made as needed

more pressing though is actual sample procurement - I
need fastqs for

```
2343 (Flowers)
2344 (Flowers)
2932 (Jang and Ehrenreich - https://www.ncbi.nlm.nih.gov/biosample/1057854)
3071
3086
GB119 (Rory)

1691 (Gallaher)
2935 (Flowers)
3059
3062
2931 (Flowers)
2342 (Flowers)
1952 (Flowers)
```

need to get these from
- Rory's 2019 paper, ENA accession PRJEB33012
- Flowers 2015
- Jang and Ehrenreich (CC2932)
- Gallaher 2015

whereas BAMs for the 30xx samples are in
`/research/data/chlamydomonas/quebec/Individual.InDelRealigned.BAMs` -
these need to be converted to fastq, perhaps via `bedtools` or
something similar

```bash
# test
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR337/004/ERR3378074/ERR3378074_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3378074/GB119.F.fq.gz
```

the two are the same after a few Python checks, minus the accession number being
added to the read id - going to download the latter files

```bash
# GB119
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3378074/GB119.F.fq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3378074/GB119.R.fq.gz

# 1691
https://www.ncbi.nlm.nih.gov/sra/SRX872474[accn]
```

## 10/8/2021

using the SRA toolkit to get these reads:

```bash
# after symlinking fastq-dump and prefetch into bin/

./bin/prefetch -v SRR1734615 # CC2931
./bin/fastq-dump \ 
--outdir /scratch/projects/chlamydomonas/recombination/parental_fastq \
--split-files ~/ncbi/public/sra/SRR1734615.sra
# renamed to CC2931_1 and _2 and gzipped - takes 55 min! 

./bin/prefetch -v SRR1797948 # CC1691
./bin/fastq-dump \
--outdir /scratch/projects/chlamydomonas/recombination/parental_fastq/ \
--split-files ~/ncbi/public/sra/SRR1797948.sra
```

other accessions I'm repeating this for:

```
SRR1734614 # CC-2344
SRR1734613 # CC-1952
SRR1734603 # CC-2935
SRR1734601 # CC-2343
SRR1734599 # CC-2342
```

before I do that - these already exist in `/research/data/chlamydomonas/species_wide/`!
that saves a good bit of time - symlinking these instead

I just need to keep 1691 and GB119 - all the Flowers strains are otherwise in this folder -
and get the 30xx strains from the path Rob sent (`/research/data/chlamydomonas/quebec/fastq/`)

also - looks like the strange case of 2932 continues, where it turns out it was technically
two fastqs (one 50 bp and one 90 bp)

## 12/8/2021

today - modify read qc + trimming workflow for parental fastqs 

fastqc will likely be its own step, after which based on the reports
I might have to do sample specific trimming - we shall see
    
```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk fastqc
```

## 13/8/2021

looks good - now to unzip:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk unzip_fastqc
```

checking the reports:

```bash
cd data/alignments/parental_fastq/fastqc/
grep -P 'FAIL\t' */summary.txt
```

got some per base sequence quality failures - that's no good - no adapter
content failures at least!

```
  1 CC1691_1_fastqc/summary.txt:FAIL        Per base sequence quality       CC1691_1.fastq.gz
  2 CC1691_1_fastqc/summary.txt:FAIL        Per sequence GC content CC1691_1.fastq.gz
  3 CC1691_2_fastqc/summary.txt:FAIL        Per base sequence quality       CC1691_2.fastq.gz
  4 CC1691_2_fastqc/summary.txt:FAIL        Per sequence GC content CC1691_2.fastq.gz
  5 CC2342_1_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2342_1.fastq.gz
  6 CC2342_2_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2342_2.fastq.gz
  7 CC2343_2_fastqc/summary.txt:FAIL        Per base sequence quality       CC2343_2.fastq.gz
  8 CC2343_2_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2343_2.fastq.gz
  9 CC2344_2_fastqc/summary.txt:FAIL        Per base sequence quality       CC2344_2.fastq.gz
 10 CC2344_2_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2344_2.fastq.gz
 11 CC2931_1_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2931_1.fastq.gz
 12 CC2932_50_1_fastqc/summary.txt:FAIL     Per base sequence content       CC2932_50_1.fastq.gz
 13 CC2932_50_2_fastqc/summary.txt:FAIL     Per base sequence quality       CC2932_50_2.fastq.gz
 14 CC2932_50_2_fastqc/summary.txt:FAIL     Per base sequence content       CC2932_50_2.fastq.gz
 15 CC2932_90_1_fastqc/summary.txt:FAIL     Per base sequence content       CC2932_90_1.fastq.gz
 16 CC2932_90_1_fastqc/summary.txt:FAIL     Per sequence GC content CC2932_90_1.fastq.gz
 17 CC2932_90_2_fastqc/summary.txt:FAIL     Per base sequence content       CC2932_90_2.fastq.gz
 18 CC2932_90_2_fastqc/summary.txt:FAIL     Per sequence GC content CC2932_90_2.fastq.gz
 19 CC2935_1_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2935_1.fastq.gz
 20 CC2935_2_fastqc/summary.txt:FAIL        Per base sequence quality       CC2935_2.fastq.gz
 21 CC2935_2_fastqc/summary.txt:FAIL        Per tile sequence quality       CC2935_2.fastq.gz
 22 CC3086_2_fastqc/summary.txt:FAIL        Per base sequence quality       CC3086_2.fastq.gz
 23 CC3086_2_fastqc/summary.txt:FAIL        Per tile sequence quality       CC3086_2.fastq.gz
 24 GB119_2_fastqc/summary.txt:FAIL Per tile sequence quality       GB119_2.fastq.gz
```

after looking at some of these `fastqc_data.txt` files, I think I ought to do a leading/trailing
trim on top of the standard sliding window trim - the window size should also still be 4 bp, 
PHRED 20 

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk trim_reads
```

## 15/8/2021

done in just about 9.5 hours without an apparent hitch - now to write rules for
fastqc and unzipping said fastqc results

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk fastqc_trim
```

## 16/8/2021

unzipping the fastqc output:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk unzip_fastqc_trim
```

checking the reports:

```bash
cd data/alignments/parental_fastq_trim/fastqc/
grep -P 'FAIL\t' */summary.txt
```

looks like there are still some standouts:

```
CC2342_trim_1_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2342_trim_1.fq.gz
CC2342_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2342_trim_2.fq.gz
CC2343_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2343_trim_2.fq.gz
CC2344_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2344_trim_2.fq.gz
CC2931_trim_1_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2931_trim_1.fq.gz
CC2932_50_trim_1_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_50_trim_1.fq.gz
CC2932_50_trim_2_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_50_trim_2.fq.gz
CC2932_90_trim_1_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_90_trim_1.fq.gz
CC2932_90_trim_2_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_90_trim_2.fq.gz
CC2935_trim_1_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2935_trim_1.fq.gz
CC2935_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2935_trim_2.fq.gz
CC3086_trim_2_fastqc/summary.txt:FAIL   Per base sequence quality       CC3086_trim_2.fq.gz
CC3086_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC3086_trim_2.fq.gz
GB119_trim_2_fastqc/summary.txt:FAIL    Per tile sequence quality       GB119_trim_2.fq.gz

```

manually looking at the per tile plots in the html files, it seems it's no more than a 
single tile each time - going to proceed as normal with these

the per base sequence content issues seem to step from a huge G-T disparity in the first
base for all four of those examples

that 3086 sample on the other hand looks like a total trainwreck - I'm going to redo
trimming with a sliding window of 30 just for that, honestly - plus I think I'd rather
this dataset have fewer but higher quality SNPs to begin with

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk trim_reads
# SLIDINGWINDOW changed to 5:30
```

8 hours later - fastqc again:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk fastqc_trim
```

## 17/8/2021

unzipping:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk unzip_fastqc_trim
```

you know the drill:

```bash
cd data/alignments/parental_fastq_trim/fastqc/
grep -P 'FAIL\t' */summary.txt
```

cleaned up some failures, but not all:

```
CC2343_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2343_trim_2.fq.gz
CC2344_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2344_trim_2.fq.gz
CC2931_trim_1_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2931_trim_1.fq.gz
CC2932_50_trim_1_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_50_trim_1.fq.gz
CC2932_50_trim_2_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_50_trim_2.fq.gz
CC2932_90_trim_1_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_90_trim_1.fq.gz
CC2932_90_trim_2_fastqc/summary.txt:FAIL        Per base sequence content       CC2932_90_trim_2.fq.gz
CC2935_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC2935_trim_2.fq.gz
CC3086_trim_2_fastqc/summary.txt:FAIL   Per base sequence quality       CC3086_trim_2.fq.gz
CC3086_trim_2_fastqc/summary.txt:FAIL   Per tile sequence quality       CC3086_trim_2.fq.gz
```

looking at the plots in the html file and 3086 still seems lower qual than I'd like
the 3' of the read despite the more aggressive trimming - I'm going to go ahead with
SNP calls and see if the GQ scores bear this out
