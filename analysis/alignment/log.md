
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

## 7/9/2021

hello from Vancouver! let's get this parental alignment/variant calling train rolling - 
need to generate bams this time around

```bash
mkdir -p data/alignments/parental_bam_temp/
mkdir -p data/alignments/parental_bam/

time snakemake -pr --cores 16 \
-s analysis/alignment/parental_alignment.smk
```

## 8/9/2021

done in 16 hours - nice! forgot to index the bams though:

```bash
time snakemake -pr -s analysis/alignment/parental_alignment.smk bam_idx
```

now for variant calling - will need to create reference VCFs for each
pair of parents - will do this in `analysis/genotyping/log.md`

## 14/9/2021

why are the trimmed 1691 files empty? time for some detective work

it seems the original fastqs are 2.7G each - but *all* reads failed
given the log:

```
TrimmomaticPE: Started with arguments:
 -threads 1 -phred33 data/alignments/parental_fastq/CC1691_1.fastq.gz data/alignments/parental_fastq/CC1691_2.fastq.gz data/alignments/parental_fastq_trim/CC1691_trim_1.fq.gz data/alignments/parental_fastq_trim/CC1691_trim_unpaired_1.fq.gz data/alignments/parental_fastq_trim/CC1691_trim_2.fq.gz data/alignments/parental_fastq_trim/CC1691_trim_unpaired_2.fq.gz SLIDINGWINDOW:5:30 LEADING:5 TRAILING:5
Input Read Pairs: 47423304 Both Surviving: 0 (0.00%) Forward Only Surviving: 0 (0.00%) Reverse Only Surviving: 0 (0.00%) Dropped: 47423304 (100.00%)
TrimmomaticPE: Completed successfully
```


let's try running trimmomatic specifically on these with less stringent parameters
(4:20 instead of 5:30 for the sliding window, for starters)

```bash
mkdir -p 1691_test

time trimmomatic PE -threads 16 -phred33 \
data/alignments/parental_fastq/CC1691_1.fastq.gz \
data/alignments/parental_fastq/CC1691_2.fastq.gz \
1691_test/CC1691_trim_1.fq.gz \
1691_test/CC1691_trim_unpaired_1.fq.gz \
1691_test/CC1691_trim_2.fq.gz \
1691_test/CC1691_trim_unpaired_2.fq.gz \
SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5  # took 15 min
```

seems most reads survived - 

```
Input Read Pairs: 47423304 
Both Surviving: 46061202 (97.13%) 
Forward Only Surviving: 1302018 (2.75%) 
Reverse Only Surviving: 46449 (0.10%) 
Dropped: 13635 (0.03%)
```
let's try this again and incrementally increase the values - 

```bash
time trimmomatic PE -threads 16 -phred33 \
data/alignments/parental_fastq/CC1691_1.fastq.gz \
data/alignments/parental_fastq/CC1691_2.fastq.gz \
1691_test/CC1691_trim_1.fq.gz \
1691_test/CC1691_trim_unpaired_1.fq.gz \
1691_test/CC1691_trim_2.fq.gz \
1691_test/CC1691_trim_unpaired_2.fq.gz \
SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 
```

also looks good -

```
Input Read Pairs: 47423304 
Both Surviving: 46067989 (97.14%) 
Forward Only Surviving: 1307487 (2.76%) 
Reverse Only Surviving: 36269 (0.08%) 
Dropped: 11559 (0.02%)
TrimmomaticPE: Completed successfully
```

## 15/9/2021

5:25?

```
Input Read Pairs: 47423304 
Both Surviving: 23471503 (49.49%) 
Forward Only Surviving: 13761704 (29.02%) 
Reverse Only Surviving: 4820153 (10.16%) 
Dropped: 5369944 (11.32%)
TrimmomaticPE: Completed successfully # done in 4 min
```

that's a pretty substantial drop off - went from 97% to 49%! I can see how
5:30 nuked it entirely

checking the ms again, and I can see why - these reads are 2x50 bp! let's
stick with the 5:25 dataset then - fewer reads but the reads need to be
the highest quality possible

replacing the trim fastqs in `data/alignment/parental_fastq_trim` with
these and rerunning the alignment workflow

```bash
rm -v data/alignments/parental_fastq_trim/CC1691_trim*
mv -v 1691_test/CC1691_* data/alignments/parental_fastq_trim/
rmdir 1691_test

time snakemake -pr --cores 16 \
-s analysis/alignment/parental_alignment.smk
```

during the dry run, snakemake detected that it needed to remake `CC1691.bam`
automatically because the 'input files [were] updated by another job' - neat! 

that said the alignment itself seems to be breaking... it keeps hitting a read
with a corrupted QUAL string, and finds two reads that are paired but have
different names for some reason - looks like trimmomatic has introduced something
funny here

```
[mem_sam_pe] paired reads have different names: "SRR1797948.47307894", "SRR1797948.47306786"

Error in job bwa_aln while creating output file data/alignments/parental_bam_temp/CC1691.bam.
RuleException:
CalledProcessError in line 39 of /research/projects/chlamydomonas/genomewide_recombination/rcmb/analysis/alignment/parental_alignment.smk:
Command 'time bwa mem -t 16 data/references/CC4532.w_organelles_MTplus.fa data/alignments/parental_fastq_trim/CC1691_trim_1.fq.gz data/alignments/parental_fastq_trim/CC1691_trim_2.fq.gz | samtools view -Sb - > data/alignments/parental_bam_temp/CC1691.bam' returned non-zero exit status 1.
  File "/research/projects/chlamydomonas/genomewide_recombination/rcmb/analysis/alignment/parental_alignment.smk", line 39, in __rule_bwa_aln
  File "/home/hasans11/.conda/env/work/lib/python3.6/concurrent/futures/thread.py", line 56, in run
```

going to try switchng the trim to 5:20 and using those files

```
Input Read Pairs: 47423304 
Both Surviving: 46067989 (97.14%) 
Forward Only Surviving: 1307487 (2.76%) 
Reverse Only Surviving: 36269 (0.08%) 
Dropped: 11559 (0.02%)
```

still breaking! looks like this problematic read is probably in the original:

```
[mem_sam_pe] paired reads have different names: "SRR1797948.47306876", "SRR1797948.47306774"

[E::sam_parse1] CIGAR and query sequence are of different length
[W::sam_read1] Parse error at line 90567849
```

going to try to generate a sam with the original - if it doesn't break, then
it's likely trimmomatic is the culprit here

```bash
time bwa mem -t 16 data/references/CC4532.w_organelles_MTplus.fa \
data/alignments/parental_fastq/CC1691_1.fastq.gz \
data/alignments/parental_fastq/CC1691_2.fastq.gz > test.sam

## | samtools view -Sb - > data/alignments/parental_bam_temp/CC1691.bam
```

## 16/9/2021

looks like this worked - those weirdly named pair names are still a thing but
the alignment completes successfully in ~20 min and the sam -> bam conversion
happens without a hitch in another 10

```
[mem_sam_pe] paired reads have different names: "SRR1797948.47308282", "SRR1797948.47306792"
```

going to try generating 4:20 reads in this case, and then trying to align with
those since it seems that the trimmed reads are experiencing some disparity
between reads and equivalent cigar strings

```bash
time trimmomatic PE -threads 16 -phred33 \
data/alignments/parental_fastq/CC1691_1.fastq.gz \
data/alignments/parental_fastq/CC1691_2.fastq.gz \
1691_test/CC1691_trim_1.fq.gz \
1691_test/CC1691_trim_unpaired_1.fq.gz \
1691_test/CC1691_trim_2.fq.gz \
1691_test/CC1691_trim_unpaired_2.fq.gz \
SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 
```

alignment time:

```bash
time bwa mem -t 20 data/references/CC4532.w_organelles_MTplus.fa \
1691_test/CC1691_trim_1.fq.gz \
1691_test/CC1691_trim_2.fq.gz > test_trim.sam # 16 min

# looks good - moment of truth
samtools view -Sb test_trim.sam > test_trim.bam
```

caused a parse error again! I don't get it

```
[W::sam_read1] Parse error at line 91308249
[main_samview] truncated file.
```

looking at this line directly in python - it seems the read name is lopped off!

```
SRR1797948.46814928     163     cpDNA   202184  60      101M    =       202839  756     TAGTTACCCGAAGGGGTTTACATACTCCGAAGGAGGAAGCAGGCAGTGGTACAATAAATAACAATAAATAACAATAAATAACTAGTATATAAATATAGGAT     ;9:999;<<;9:<<;;:889;999:::;;::<<:<<99<<:<;<:<:<<:9:::99::99:;::99::99:;;:9:::9::;;;=;:9989:;9999<:89     NM:i:0  MD:Z:101 MC:Z:101M        AS:i:101        XS:i:47

TCCGGGGAGTGCGTTGACTCCAGAAGGTTATTTGCCCGTCCACGATTG        ;:9;99<:<;;9:<;;;;9;:;<;:9;:;::;;<;9<<9999::<<;;;;;<:;;999<       NM:i:0  MD:Z:59 MC:Z:101M       AS:i:59 XS:i:38
```

let's try removing this and then converting to bam:

```python
from tqdm import tqdm
with open('test_fixed.sam', 'w') as f_out:
    with open('test_trim.sam', 'r') as f_in:
        counter = 0
        for line in tqdm(f_in):
            if line.startswith('@') or line.startswith('SRR'):
                f_out.write(line)
            else:
                counter += 1
```

converting:

```bash
samtools view -Sb test_fixed.sam > test_fixed.bam
```

this worked! hooray! going to have to update snakemake with
some sort of stipulation for this - either that or remove 1691
from it entirely, and keep a separate workflow for it

```bash
# trying with the 5:25 trimmed reads

time bwa mem -t 20 data/references/CC4532.w_organelles_MTplus.fa \
data/alignments/parental_fastq_trim/CC1691_trim_1.fq.gz \
data/alignments/parental_fastq_trim/CC1691_trim_2.fq.gz > CC1691.sam # 16 min
```

doesn't seem reads are bonked the same way:

```python
>>> from tqdm import tqdm
>>> with open('CC1691_fixed.sam', 'w') as f_out:
...     counter = 0
...     with open('CC1691.sam', 'r') as f:
...         for line in tqdm(f):
...             if line.startswith('@') or line.startswith('SRR'):
...                 f_out.write(line)
...             else:
...                 counter += 1
...                 continue
91724748it [02:29, 612131.93it/s]
>>> counter
0
```

let's try the conversion:

```bash
time samtools view -Sb CC1691.sam > CC1691.bam
```

this worked?? let's try snakemake again then:

```bash
time snakemake -pr -s analysis/alignment/parental_alignment.smk \
--cores 20
```

## 17/9/2021

it broke again somehow - I'm so confused

at this point I'm going to split the sam generation and sam -> bam
conversion into two steps


```bash
time snakemake -pr -s analysis/alignment/parental_alignment.smk \
--cores 20
# now with new rule bam_convert
```

so this breaks - but when I run the command independently in the shell
it's fine? I'll just do that and have snakemake pick up from `bam_convert`
in this case

```bash
time bwa mem -t 20 data/references/CC4532.w_organelles_MTplus.fa \
data/alignments/parental_fastq_trim/CC1691_trim_1.fq.gz \
data/alignments/parental_fastq_trim/CC1691_trim_2.fq.gz > \
data/alignments/parental_bam_temp/CC1691.sam

time snakemake -pr --cores 20 -s analysis/alignment/parental_alignment.smk
```

THIS broke??? after it worked yesterday independent of snakemake?

let's scan the outfile:

```python
from tqdm import tqdm
with open('data/alignments/parental_bam_temp/CC1691.sam', 'r') as f:
    counter = 0
    lines = []
    for i, line in tqdm(enumerate(f)):
        if line.startswith('SRR') or line.startswith('@'):
            continue
        else:
            counter += 1
            lines.append(i)
            lines.append(line)
```

looks like this is the culprit:

```
>>> lines
[91643941, '9<799;:;;;=:;:9<<;<<=;<>:>:9<:<=<=99<<;<=;<;<;<;==>==::;;>>=>:9:;\tAS:i:0\tXS:i:0\n']
```

so let's create a corrected file without that truncated read and get
back to business:

```python
from tqdm import tqdm
with open('data/alignments/parental_bam_temp/CC1691_fixed.sam', 'w') as f_out:
    with open('data/alignments/parental_bam_temp/CC1691.sam', 'r') as f:
        counter = 0
        lines = []
        for i, line in tqdm(enumerate(f)):
            if line.startswith('SRR') or line.startswith('@'):
                f_out.write(line)
            else:
                counter += 1
                lines.append(i)
                lines.append(line)
```

do some renaming and get snakemake going again:

```bash
mv -v data/alignments/parental_bam_temp/CC1691.sam \
data/alignments/parental_bam_temp/CC1691_old.sam

mv -v data/alignments/parental_bam_temp/CC1691_fixed.sam \
data/alignments/parental_bam_temp/CC1691.sam

time snakemake -pr --cores 20 -s analysis/alignment/parental_alignment.smk
```

it worked! done in 2 hours - I ought to create a 'sam repair' script 
for reproducibility's sake down the line

## 7/10/2021

and we're back - now to fix the CC2932 issue, where the
two bams need to be combined

updated the PREFIXES variable to remove the `_50` and `_90` prefixes,
and added a `combine_2932` rule that specifically concatenates those two
into `CC2932.sam`

```bash
time snakemake -pr -s analysis/alignment/parental_alignment.smk --cores 20
```

## 20/10/2021

finally - time to prep bams for phase change analysis using readcomb! 

testing on a single bam:

```bash
mkdir data/alignments/bam_filtered

time readcomb-bamprep --bam data/alignments/bam/2343x1691.bam \
--threads 16 --outdir data/alignments/bam_filtered/
```

took a bit of debugging but looks good - now to add this to the
alignment workflow

```bash
# updated with readcomb_bamprep rule
time snakemake -pr -s analysis/alignment/alignment.smk --cores 16
```

## 23/10/2021

took 5 hours but files look good

getting `readcomb-filter` going -

```bash
mkdir -p data/phase_changes/sam/
mkdir -p data/phase_changes/bam/

time readcomb-filter \
--bam data/alignments/bam_filtered/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--processes 8 \
--log filter_test.log \
--quality 30 \
--out data/phase_changes/bam/2343x1691.filtered.bam
```

## 25/10/2021

this bugged out due to a logging error - trying again -
 
```bash
time readcomb-filter \
--bam data/alignments/bam_filtered/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--processes 8 \
--log filter_test.log \
--quality 30 \
--out data/phase_changes/sam/2343x1691.filtered.sam
```

## 26/10/2021

looks good! looks like 13% of reads 'passed':

```
3565309 phase change-containing read pairs from total 26238981 read pairs
8827201 reads had no-match variants
8092147 reads did not have enough variants (> 0) to call
```

and this completed in 110 min given ten processes - which bodes really well

time to get this workflow going - going to continue this in `analysis/phase_changes/log.md`
(which is where I should have been in the first place)

## 7/1/2022

back here cause `CC1691.bam` continues to be a total jerk - this time since there
seems to be a read in chr3 that passes all flags but has no mate (I wonder if the mate
is that corrupted read I removed earlier?) 

either way, I'm going to have to do some manual fixes to toss this read into the gutter as well - 
from earlier, but updated:

```python
from tqdm import tqdm
with open('data/alignments/parental_bam_temp/CC1691_fixed.sam', 'w') as f_out:
    with open('data/alignments/parental_bam_temp/CC1691_old.sam', 'r') as f:
        trash_read = 'SRR1797948.46982639'
        counter = 0
        lines = []
        for i, line in tqdm(enumerate(f)):
            if line.startswith('SRR') or line.startswith('@') and not trash_read in line:
                f_out.write(line)
            else:
                counter += 1
                lines.append(i)
                lines.append(line)
```

once more, with feeling:

```bash
mv -v data/alignments/parental_bam_temp/CC1691.sam \
data/alignments/parental_bam_temp/CC1691_old.sam

mv -v data/alignments/parental_bam_temp/CC1691_fixed.sam \
data/alignments/parental_bam_temp/CC1691.sam

# will redo workflow specifically for CC1691
time snakemake -pr --cores 20 -s analysis/alignment/parental_alignment.smk
```

and then back to the `phase_changes` log we go to generate more false positives

aaaand update: the read is _still_ there - just going to remove it manually at this point,
I don't have the goddamn brainpower to reverse engineer how it managed to stick around
despite me removing it earlier

```
import pysam
from tqdm import tqdm

reader = pysam.AlignmentFile('data/alignments/parental_bam_filtered/CC1691.sorted.bam', 'rb')
writer = pysam.AlignmentFile(
    'data/alignments/parental_bam_filtered/CC1691.fixed.bam', 'wh', template=reader)

total, counter = 0, 0
for record in tqdm(reader):
    total += 1
    if record.reference_name == 'chromosome_03' and record.query_name == 'SRR1797948.46982639':
        continue
    else:
        counter += 1
        writer.write(record)
# >>> counter
# 82377392
# >>> total
# 82377393
```

will need to sort, recreate the index file, and then sort by read name:

```bash
# in data/alignments/parental_bam_filtered
rm CC1691.sorted.bai # remove old index

# sort for index generation
time samtools sort -O bam -@12 -o CC1691.fixed.sorted.bam --verbosity 4 CC1691.fixed.bam

# regen index
time samtools index CC1691.fixed.sorted.bam CC1691.sorted.bai

# create read name sorted file
rm CC1691.sorted.bam # delete old version
time samtools sort -n -O bam -@12 -o CC1691.sorted.bam CC1691.fixed.sorted.bam
```

## 10/1/2022

getting some coverage stats for the methods section

```python
counts = [40066453, 26544459, 30459259, 22726219, 33717850, 51466264, 
    22966295, 27814662, 19991410, 29707153, 30581007, 24920406, 22568234, 
    20106348, 20172916, 26150662, 20602843, 27399502, 21516400, 19754383, 
    18106544, 48260547, 20410697, 27177597, 33705747, 43510346, 20323566, 
    22981901, 22980268, 28922949] # from the excel sheet

sum(counts) / len(counts) # ~27m
```

## 20/1/2022

today - trying out the HC -bamout argument to see if I can resolve
sketchy paralogous regions

let's try this with 2932x3062 since that's been the guinea pig cross for
a bit now

going to run this on CC2932 alone to start with - otherwise maintaining the same parameters

```bash
mkdir -p data/hc_test

time gatk --java-options "-Xmx4g" HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC2932.bam \
-O data/hc_test/CC2932.vcf.gz \
-bamout data/hc_test/CC2932.bamout.bam \
-L chromosome_01 -ploidy 2 --heterozygosity 0.02 \
--indel-heterozygosity 0.002
# took 1 hr
```

I haven't set the -ERC to GVCF since I don't really care about the VCF here

seems the bamout file only has 2036779 'non HC' reads, far fewer than the original bam -
the reconstructed haplotypes are kept as 'HC' reads

making an HC-read only version for IGV:

```python
import pysam
from tqdm import tqdm

reader = pysam.AlignmentFile('data/hc_test/CC2932.bamout.bam', 'rb')
writer = pysam.AlignmentFile('data/hc_test/CC2932.hc_only.sam', 'wh', template=reader)

for record in tqdm(reader):
    if record.query_name.startswith('HC'):
        writer.write(record)

# and then compress and samtools index the output sam
```

## 12/4/2022

the day has finally come - 2 x 250 bp for all the parents! 

going to clear out the old parental bams, since I don't need these anymore.
for now - 

```bash
# in data/alignments
mkdir parental_old
mv -v parental_* parental_old/
```

and now to download the files - going to have details on Notion

```bash
# downloaded to /scratch/projects/chlamydomonas/recombination/other_stuff
# in that dir:
mv -v *xCC* ../fastq # two recombinant samples - 3071x2931 and 3071x3062
mv -v *CC* ../parental_fastq # the rest of the samples 
mv -v *GB119* ../parental_fastq # and the one GB sample
```

comparing checksums - need to download to local machine and then transfer to server

```bash
# in parental_fastq
time md5sum -c readSets.md5
```

this will raise some issues given that I've moved some of the files over to `other_stuff` - 
will check those separately

current samples look good -

```
NS.1843.002.IDT_i7_13---IDT_i5_13.CC2931_R1.fastq.gz: OK
NS.1843.002.IDT_i7_13---IDT_i5_13.CC2931_R2.fastq.gz: OK
NS.1843.002.IDT_i7_1---IDT_i5_1.CC2344_R1.fastq.gz: OK
NS.1843.002.IDT_i7_1---IDT_i5_1.CC2344_R2.fastq.gz: OK
NS.1843.002.IDT_i7_181---IDT_i5_181.CC2935_R1.fastq.gz: OK
NS.1843.002.IDT_i7_181---IDT_i5_181.CC2935_R2.fastq.gz: OK
NS.1843.002.IDT_i7_109---IDT_i5_109.CC1952_R1.fastq.gz: OK
NS.1843.002.IDT_i7_109---IDT_i5_109.CC1952_R2.fastq.gz: OK
NS.1843.002.IDT_i7_98---IDT_i5_98.CC3059_R1.fastq.gz: OK
NS.1843.002.IDT_i7_98---IDT_i5_98.CC3059_R2.fastq.gz: OK
NS.1843.002.IDT_i7_133---IDT_i5_133.CC2343_R1.fastq.gz: OK
NS.1843.002.IDT_i7_133---IDT_i5_133.CC2343_R2.fastq.gz: OK
NS.1843.002.IDT_i7_121---IDT_i5_121.CC2342_R1.fastq.gz: OK
NS.1843.002.IDT_i7_121---IDT_i5_121.CC2342_R2.fastq.gz: OK
NS.1843.002.IDT_i7_169---IDT_i5_169.CC2932_R1.fastq.gz: OK
NS.1843.002.IDT_i7_169---IDT_i5_169.CC2932_R2.fastq.gz: OK
NS.1843.002.IDT_i7_97---IDT_i5_97.CC1691_R1.fastq.gz: OK
NS.1843.002.IDT_i7_97---IDT_i5_97.CC1691_R2.fastq.gz: OK
NS.1843.002.IDT_i7_122---IDT_i5_122.CC3071_R1.fastq.gz: OK
NS.1843.002.IDT_i7_122---IDT_i5_122.CC3071_R2.fastq.gz: OK
NS.1843.002.IDT_i7_134---IDT_i5_134.CC3086_R1.fastq.gz: OK
NS.1843.002.IDT_i7_134---IDT_i5_134.CC3086_R2.fastq.gz: OK
NS.1843.002.IDT_i7_110---IDT_i5_110.CC3062_R1.fastq.gz: OK
NS.1843.002.IDT_i7_110---IDT_i5_110.CC3062_R2.fastq.gz: OK
NS.1843.002.IDT_i7_146---IDT_i5_146.GB119_R1.fastq.gz: OK
NS.1843.002.IDT_i7_146---IDT_i5_146.GB119_R2.fastq.gz: OK
```

going to make a copy in `../fastq` just for those two recombinant samples that have
since been added:

```bash
cp -v readSums.md5 ../fastq
cd ../fastq
time md5sum -c readSets.md5
```

also looks good! 

```
NS.1843.002.IDT_i7_158---IDT_i5_158.CC3071xCC2931_R1.fastq.gz: OK
NS.1843.002.IDT_i7_158---IDT_i5_158.CC3071xCC2931_R2.fastq.gz: OK
NS.1843.002.IDT_i7_170---IDT_i5_170.CC3071xCC3062_R1.fastq.gz: OK
NS.1843.002.IDT_i7_170---IDT_i5_170.CC3071xCC3062_R2.fastq.gz: OK
```

time to symlink and then rename said symlinks to something simpler

```bash
# in data/alignments/fastq
ln -sv /scratch/projects/chlamydomonas/recombination/fastq/NS.1843* .

# then renamed symlinks to follow prev format - not going to do python, just four files
```

QC and trimming is next - this should be straightforward for the parents, but
not the recombinants. `read_qc_trim.smk` wasn't written with outputs explicitly
specified, or even a rule all - I'm just going to keep this simple and run the
commands directly with a temp bash script for the recombinants, but otherwise
use the existing parental read qc snakefile (with minor adjustments, since the other
had some edge case corrections) for the parents

let's also symlink the parents while we're at it

```bash
# in data/alignments/parental_fastq
ln -sv /scratch/projects/chlamydomonas/recombination/parental_fastq/NS*fastq.gz .
```

and renaming these in Python, as before - o

```python
import os
import re
import glob

fname_list = glob.glob('*R?.fastq.gz')

pattern = '[CCGB]{2}[0-9]{3,4}_R[12].fastq.gz'
for fname in fname_list:
    new_fname = re.search(pattern, fname).group()
    os.rename(fname, new_fname)
```

tomorrow - run a bash script to QC the two new recombinant sequences, and then
update + run the parental fastq QC workflow as well

## 13/4/2022

QCing the recombinant sequences first - 

```bash
# take 25 min each
time fastqc --kmers 7 --outdir data/alignments/fastq/fastqc/raw \
data/alignments/fastq/3071x2931_R*.fastq.gz

time fastqc --kmers 7 --outdir data/alignments/fastq/fastqc/raw \
data/alignments/fastq/3071x3062_R*.fastq.gz
```

unzipping the fastqc files:

```bash
unzip data/alignments/fastq/fastqc/raw/3071x2931_R1_fastqc.zip \
-d data/alignments/fastq/fastqc
unzip data/alignments/fastq/fastqc/raw/3071x2931_R2_fastqc.zip \
-d data/alignments/fastq/fastqc

unzip data/alignments/fastq/fastqc/raw/3071x3062_R1_fastqc.zip \
-d data/alignments/fastq/fastqc
unzip data/alignments/fastq/fastqc/raw/3071x3062_R2_fastqc.zip \
-d data/alignments/fastq/fastqc
```

checking for failed tests:

```bash
cd data/alignments/fastq/fastqc
grep -P 'FAIL\t' 3071*/summary.txt
```

only GC content + adapter content issues - let's get a move on then and trim

```bash
for prefix in 3071x3062 3071x2931; do
    time trimmomatic PE -threads 16 -phred33 \
    data/alignments/fastq/${prefix}_R1.fastq.gz \
    data/alignments/fastq/${prefix}_R2.fastq.gz \
    data/alignments/fastq_trim/${prefix}_trim_1.fq.gz \
    data/alignments/fastq_trim/${prefix}_trim_unpaired_1.fq.gz \
    data/alignments/fastq_trim/${prefix}_trim_2.fq.gz \
    data/alignments/fastq_trim/${prefix}_trim_unpaired_2.fq.gz \
    ILLUMINACLIP:bin/NEBNext_dual.fasta:2:30:10 \
    SLIDINGWINDOW:4:20
done
```
    
## 14/4/2022

fastqc'ing the trimmed reads - 
    
```bash
time fastqc --kmers 7 --outdir data/alignments/fastq_trim/fastqc/raw \
data/alignments/fastq_trim/3071x2931_trim_*.fq.gz

time fastqc --kmers 7 --outdir data/alignments/fastq_trim/fastqc/raw \
data/alignments/fastq_trim/3071x3062_trim_*.fq.gz

# whoops - this did paired ones as well - will just remove those
```

unzipping - 

```bash
unzip data/alignments/fastq_trim/fastqc/raw/3071x2931_trim_1_fastqc.zip \
-d data/alignments/fastq_trim/fastqc
unzip data/alignments/fastq_trim/fastqc/raw/3071x2931_trim_2_fastqc.zip \
-d data/alignments/fastq_trim/fastqc

unzip data/alignments/fastq_trim/fastqc/raw/3071x2931_trim_1_fastqc.zip \
-d data/alignments/fastq_trim/fastqc
unzip data/alignments/fastq_trim/fastqc/raw/3071x2931_trim_2_fastqc.zip \
-d data/alignments/fastq_trim/fastqc
```

checking with grep command above - nothing but GC content issues in 3071 x 2931, as expected

and now for alignment! just manually updating `data/alignments/samples.txt` with the added two

```bash
time snakemake -pr -s analysis/alignment/alignment.smk --cores 16
```

while this is running: getting the parental QC workflow going

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk fastqc # took 2 hrs
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk unzip_fastqc
```

seems all the failures are GC content + adapter content - as expected

trimming:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk --cores 16 trim_reads
```

## 15/4/2022

done in 4 hrs - now to fastqc the trimmed reads:

```bash
time snakemake -pr -s analysis/alignment/parental_fastq_qc.smk fastqc_trim
```

fastqc results look good - had to add ILLUMINACLIP to the trimmomatic command but otherwise
things were in the clear

time for alignment! 

```bash
time snakemake -pr -s analysis/alignment/parental_alignment.smk --cores 20
# took 16 hours
```

back to genotyping! 

## 22/5/2022

out of curiosity - going to get some MAPQ summaries for each parental bam -

```python
>>> import pysam
>>> from tqdm import tqdm
>>> reader = pysam.AlignmentFile('data/alignments/parental_bam/CC3071.bam')
>>> scores = {}
>>> counter = 0
>>> for record in tqdm(reader):
...     if not record.mapq in scores:
...         scores[record.mapq] = 1
...     elif record.mapq in scores:
...         scores[record.mapq] += 1
...     counter += 1
...     if counter > 2e7:
...         break
19950904it [00:39, 522635.42it/s]
>>> reader = pysam.AlignmentFile('data/alignments/parental_bam/CC2342.bam')
>>> scores_2 = {}
>>> counter = 0
>>> for record in tqdm(reader):
...     if not record.mapq in scores_2:
...         scores_2[record.mapq] = 1
...     elif record.mapq in scores_2:
...         scores_2[record.mapq] += 1
...     counter += 1
...     if counter > 2e7:
...         break
18781262it [00:37, 496688.30it/s]
>>> for score in sorted(scores.keys()):
...     print(score, scores[score], scores_2[score])
0 2816164 3863909
1 33472 48701
2 28164 42281
3 36705 52358
4 28782 41085
5 27143 40919
6 35883 47273
7 27464 37593
8 17051 26817
9 29504 43105
10 13854 21906
11 12455 18584
12 17297 26711
13 15548 23354
14 11672 17869
15 21498 33397
16 16495 21904
17 18646 25813
18 22831 29790
19 25337 31995
20 25028 29419
21 29536 42220
22 35491 39408
23 24906 30424
24 35506 42967
25 43219 52995
26 10272 16331
27 67501 68798
28 11641 18015
29 10721 15797
30 18168 25232
31 12903 20215
32 10129 15908
33 21079 32236
34 10088 14716
35 11446 15373
36 13489 20942
37 14344 19135
38 9170 13377
39 22332 31654
40 272289 250652
41 22863 22834
42 25434 30585
43 28114 32579
44 36071 37621
45 37939 50594
46 67623 54689
47 23400 29798
48 20213 27772
49 21162 24691
50 16836 21496
51 17747 23179
52 24100 30448
53 17194 25207
54 27204 32437
55 20718 27768
56 24702 32149
57 40952 45138
58 54486 35229
59 24442 31858
60 15483578 12854012
```

it seems the majority of reads are in fact MAPQ 60 - which is curious given
some of these extremely scuffed regions I'm seeing in 3071x2931, which has
a ridiculous amount of phase changes to the point where I don't believe a single one

doesn't seem like filtering reads below MAPQ 60 actually has a huge impact
on things:

```bash
$ head -n 1 3071.filt.cov.txt
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq

hasans11@hpcnode1 /research/projects/chlamydomonas/genomewide_recombination/rcmb/data/alignments/parental_bam (work)
$ grep 'chromosome' 3071*txt
3071.cov.txt:chromosome_01      1       8225636 2152772 7920016 96.2845 53.6853 35.9    47.6
3071.cov.txt:chromosome_02      1       8655884 2092944 8383456 96.8527 50.2713 35.9    53.8
3071.cov.txt:chromosome_03      1       9286894 2237802 8842474 95.2145 49.6341 35.9    51.8
3071.cov.txt:chromosome_04      1       4130073 972317  3801761 92.0507 44.2681 35.9    49.3
3071.cov.txt:chromosome_05      1       3682160 895837  3341851 90.7579 46.5652 35.9    47.5
3071.cov.txt:chromosome_06      1       8913359 2028879 8219514 92.2157 46.3422 35.9    52
3071.cov.txt:chromosome_07      1       6492107 1602165 6240800 96.129  50.211  35.9    51.7
3071.cov.txt:chromosome_08      1       4526983 1252945 4207986 92.9534 55.1016 35.9    42.3
3071.cov.txt:chromosome_09      1       6807148 1561383 6229905 91.52   45.3559 35.9    49.9
3071.cov.txt:chromosome_10      1       6800247 1633219 6488750 95.4193 49.4355 35.9    51.9
3071.cov.txt:chromosome_11      1       4479522 1021650 4039110 90.1683 43.4616 35.9    50.6
3071.cov.txt:chromosome_12      1       9952739 2377667 9469546 95.1451 49.8073 35.9    51.6
3071.cov.txt:chromosome_13      1       5281438 1281443 4944631 93.6228 47.1286 35.9    49.7
3071.cov.txt:chromosome_14      1       4217303 1337111 3912525 92.7732 64.2766 35.9    37
3071.cov.txt:chromosome_15      1       5870643 964759  3844905 65.4938 30.4157 35.9    36.1
3071.cov.txt:chromosome_16      1       8042475 1936949 7648358 95.0996 49.0719 35.9    51.6
3071.cov.txt:chromosome_17      1       6954842 1613978 6513782 93.6582 46.7952 35.9    51.8
3071.filt.cov.txt:chromosome_01 1       8225636 1604084 7111782 86.4587 42.0856 35.9    60
3071.filt.cov.txt:chromosome_02 1       8655884 1809139 7841319 90.5895 45.5451 35.9    60
3071.filt.cov.txt:chromosome_03 1       9286894 1843850 8035583 86.5261 42.9557 35.9    60
3071.filt.cov.txt:chromosome_04 1       4130073 742668  3268206 79.1319 36.4288 35.9    60
3071.filt.cov.txt:chromosome_05 1       3682160 644095  2879544 78.2026 36.0736 35.9    60
3071.filt.cov.txt:chromosome_06 1       8913359 1666484 7337316 82.3182 40.088  35.9    60
3071.filt.cov.txt:chromosome_07 1       6492107 1303905 5739641 88.4095 43.0413 35.9    60
3071.filt.cov.txt:chromosome_08 1       4526983 821377  3740530 82.6274 37.4909 35.9    60
3071.filt.cov.txt:chromosome_09 1       6807148 1216424 5435904 79.8558 37.6205 35.9    60
3071.filt.cov.txt:chromosome_10 1       6800247 1348341 5914958 86.9815 42.8646 35.9    60
3071.filt.cov.txt:chromosome_11 1       4479522 802308  3542483 79.0817 36.6408 35.9    60
3071.filt.cov.txt:chromosome_12 1       9952739 1938320 8460427 85.006  42.4393 35.9    60
3071.filt.cov.txt:chromosome_13 1       5281438 1005289 4359335 82.5407 39.8962 35.9    60
3071.filt.cov.txt:chromosome_14 1       4217303 756509  3443964 81.6627 37.0256 35.9    60
3071.filt.cov.txt:chromosome_15 1       5870643 471082  2187177 37.2562 15.9642 35.9    60
3071.filt.cov.txt:chromosome_16 1       8042475 1573595 6879314 85.5373 42.1232 35.9    60
3071.filt.cov.txt:chromosome_17 1       6954842 1313633 5822915 83.7246 40.2566 35.9    60
```

as a test, I'm going to create filtered versions of 3071 and 2931 and specifically redo
the phase change filtering for that cross

going to also do the same for 2344 and 2931 to see if that has any impact - since 2344x2931
seems to have a reasonable set of phase changes while 3071x2931 is straight from the bowels
of hell currently

will do these tests in a new temp dir called `mapq_test`

```bash
mkdir mapq_test

samtools view -O bam -q 60 data/alignments/parental_bam/CC3071.bam > mapq_test/CC3071.60.bam
samtools view -O bam -q 60 data/alignments/parental_bam/CC2931.bam > mapq_test/CC2931.60.bam
samtools view -O bam -q 60 data/alignments/parental_bam/CC2344.bam > mapq_test/CC2344.60.bam

# in mapq_test
for fname in *bam; do samtools index ${fname}; done

# this should technically be in genotyping/ but what the hell
parallel -j5 'time freebayes -f data/references/CC4532.w_organelles_MTplus.fa \
--theta 0.02 --ploidy 2 --genotype-qualities -r chromosome_{} \
mapq_test/CC3071.60.bam mapq_test/CC2931.60.bam > mapq_test/3071x2931.chr{}.vcf
' ::: {01..09} {10..17}

parallel -j2 'time freebayes -f data/references/CC4532.w_organelles_MTplus.fa \
--theta 0.02 --ploidy 2 --genotype-qualities -r {} \
mapq_test/CC3071.60.bam mapq_test/CC2931.60.bam > mapq_test/3071x2931.{}.vcf
' ::: mtDNA cpDNA
```

## 23/5/2022

ran overnight - now to concat these files


```
# in mapq_test
for fname in *vcf; do bgzip ${fname}; tabix ${fname}.gz; done
bcftools concat 3071x2931.*.vcf.gz > 3071x2931.vcf
bgzip 3071x2931.vcf
tabix 3071x2931.vcf.gz
```

hold on - I looked at the original 3071 file (from a prev study) and
it explains a lot of these weird calls - the 2 x 250 bp 3071 has a ton
of strange new variants that are causing false recombination calls

I'm going to check and see whether 3086 has a similar thing going, and if so,
going to redo variant calling with the older 3071 and 3086 before redoing
phase change detection for both as well

back to the `genotyping` log









