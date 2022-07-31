
## 3/6/2022

today - installing LDhelmet 1.10 on the server and testing it using
the provided example file

```bash
cd ~/apps
wget https://github.com/popgenmethods/LDhelmet/archive/refs/heads/main.zip
# modified Makefile to add libraries 
# /research/tmp_apps/my_boost/include and /research/tmp_apps/my_gsl/include
make
```

running the example workflow - to do once phase changes and such are done,
since I'm hammering the server enough as it stands

```bash
cd example_scripts
time bash run_all.bash
```

## 29/7/2022

finally getting back to work on this - need to generate parental fastas
to then feed as input to LDhelmet

seems the best way to do this is to use `bcftools concat` on a GVCF - though
for this I'll have to remake the 3071 GVCF, but that's about it (since the 2x250
one was totally bonked, and I had to revert to 2x150)

let's first try this on a GVCF that I know is good to go:

```bash
time bcftools consensus \
--fasta-ref data/references/CC4532.w_organelles_MTplus.fa \
--output 1691_consensus_test.fasta \
data/genotyping/gvcf_sample/CC1691.g.vcf.gz
```

this looks good to go, and runs super quick (done in 20 sec!) 

I'm going to get the 3071 GVCF going for now - after this, will have to make
FASTAs for each chromosome, with `chromosome_n.fasta` containing said
chromosome from each sample - should be a quick bit of SeqIO wrangling

here goes the 3071 GVCF - this will run overnight

```bash
time gatk HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC3071.bam \
-ERC GVCF -ploidy 2 \
--heterozygosity 0.02 --indel-heterozygosity 0.002 \
--native-pair-hmm-threads 2 \
-O data/genotyping/gvcf_sample/CC3071.g.vcf.gz
# took 14.5 hours
```

## 30/7/2022

today - creating fastas for each of the samples,
and then converting those fastas into per-chr
files

let's get our directories in order and then get things going:

```bash
mkdir -p data/ldhelmet/
mkdir -p data/ldhelmet/fasta
mkdir -p data/ldhelmet/fasta_chrom
```

I just realized though that I probably should
use a snakemake workflow for this whole thing,
and that includes generating the fastas

will need an ancillary script for generating the per-chr
files from the per-sample files but otherwise I should be good

```bash
# contains a rule called create_sample_fasta
# rule all currently set to output of that rule
time snakemake -pr -s analysis/ldhelmet/ldhelmet.smk # done in 7 min
```

and now to convert these fastas into per-chr files -
for this, I think I should use snakemake's script feature
and just create a small custom script for this

will need to provide all filenames at once using `expand()`
for this I think, and then update rule all with per-chr files

```bash
# updated rule all to include per-chr fastas in new dir
time snakemake -pr -s analysis/ldhelmet/ldhelmet.smk # took 10 min
```

looks great! tomorrow - set up LDhelmet rules (maybe with separate scripts?
figure this out) for NA1-only, NA2-only, and NA1xNA2 (ie all samples)



