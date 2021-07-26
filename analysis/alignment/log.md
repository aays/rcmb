
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

trimming finally done - took a day and 6 hours, annd that's with
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
time snakemake -pr -s analysis/alignment/alignment.smk --cores 4
```

fix mate info:

```bash
time snakemake -pr -s analysis/alignment/alignment.smk
```






    

