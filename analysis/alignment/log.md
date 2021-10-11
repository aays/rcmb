
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


