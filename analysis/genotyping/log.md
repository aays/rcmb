
## 8/9/2021

today - writing HaplotypeCaller snakefile for parental VCFs

```bash
touch analysis/genotyping/variant_calling.smk
```

going to be using the GATK best practices approach - 

1. individual GVCFs using HC in GVCF mode 
2. combining the GVCFs
3. joint genotyping 

going to start with 2343x1691 before generalizing using a snakefile just to make sure
this works as intended

```bash
mkdir -p data/genotyping/gvcf_sample
mkdir -p data/genotyping/gvcf_combined
mkdir -p data/genotyping/vcf

time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC2343.bam \
-ERC GVCF \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-ploidy 2 \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/genotyping/gvcf_sample/CC2343.g.vcf.gz
```

## 12/9/2021

that took 40 whole hours - oh my - might need to find some way
to parallelize this over regions down the line

but for now, just going to run 1691 to continue the 'test workflow'

```bash
time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC1691.bam \
-ERC GVCF \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-ploidy 2 \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/genotyping/gvcf_sample/CC1691.g.vcf.gz
```

## 14/9/2021

wait - this ended in like 4 minutes since the 1691 bam was never
correctly generated - what on earth? 

let's check on the fastqs - back to the alignment log I go

## 18/9/2021

and we're back after that sordid affair! 

running 1691:

```bash
time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC1691.bam \
-ERC GVCF \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-ploidy 2 \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/genotyping/gvcf_sample/CC1691.g.vcf.gz
```

## 24/9/2021

took 7 hours and we're in business! 

going to get the snakefile going - this ought to
run for a full week given it has 12 files to run through

I wonder if I can split it by chromosome - we've
got chromosomes 1 through 17, `contig_18` to `contig_57`,
and then cpDNA/mtDNA/`MT_plus_R`

## 25/9/2021

first up - bit the bullet and switched to GATK4, which I should have done
far sooner 

trying to parallelize across intervals - looks like if I'm to do it with intervals,
I need to set the intervals themselves as the 'wildcard', which only works if I'm
working with one VCF - one way to work around this _could_ be to create indiv
dirs for each sample and indiv interval list files for each interval, 
but this seems like a lot of temp-file-making for not a ton of gain

that said, if more cores are provided at the command line snakemake will start
processes on subsequent samples as they come up

setting threads to 2 and providing 6 cores - 
see [here](https://bioinformatics.stackexchange.com/questions/4608/increase-number-of-threads-for-gatk-4-0-haplotypecaller)
for more on the `--native-pair-hmm-threads` arg (apparently it doesn't use more than 2 no matter
what) as well as the implementation in [here](https://gitlab.com/computational-biology/ovarflow/-/blob/master/OVarFlow_src/Snakefile)
(roughly around line 643) 

either way - 

```bash
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 6
```

## 29/9/2021

all done - took 84 hours (e.g. 3.5 days) which is much better than I expected! 

next up - creating combined VCFs for the parental pairs and then running GenotypeGVCFs

it seems the GATK4 method is to use GenomicsDBImport instead of
CombineGVCFs to scale performance better when the number of samples increases,
but I'm going to have two samples at a time so I'll stick to CombineGVCFs

it does occur to me though that I have two CC2932s - one from the 50 bp
reads and another from the 90! going to have to make a separate
rule to handle these

```bash
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 4
# runs call_2932 rule
```

## 7/10/2021

so this breaks since GVCF mode can only handle one bam at a time

I think I need to combine the bams, but this'll require some modification of the 
alignment workflow - so back to that log I go

once that's all well and good, back to trying to generate the 2932 VCF after updating
prefixes and removing the `call_2932` rule

first, need to read the various parental pairs into snakemake somehow, since this
governs which GVCFs end up getting combined some way

here goes, starting off with the 2932 GVCF

```bash
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 10
```

## 8/10/2021

and now that that's done (took 16 hours...) for combining gvcfs:

```bash
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 10
# now with combine_parents rule that uses data/genotyping/samples.txt but commented
out
```

## 10/10/2021

combined gvcfs are also ready to go (7.7 hours later) and so on we go to the final step
of using `GenotypeGVCFs` to actually perform joint genotyping

```bash
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 8
# now with genotype_gvcfs rule
```

## 11/10/2021

took two hours and now it looks like we're in business! 

## 20/10/2021

revisiting this - I haven't yet run `readcomb-vcfprep` on the files - should
add that to the workflow and get that running! 


```bash
mkdir data/genotyping/vcf_filtered

# trying to filter a single file before adding to workflow
# readcomb 0.0.5
time readcomb-vcfprep \
--vcf data/genotyping/vcf/GB119x3062.vcf.gz \
--no_hets --min_GQ 30 \
--out data/genotyping/vcf_filtered/GB119x3062.vcf.gz
```

looks good - retained 43% of variants (1.65m from 3.79m)

updating the snakemake workflow:

```bash
# will run readcomb_vcfprep rule
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 8
```

done in half an hour

although this reminds me I need to run bamprep as well! back to the alignment log

## 2/11/2021

today - finding a quebec all samples v6 VCF to use for SNP curation - with the
idea being that if there are regions with dense het calls, these are likely
paralogous regions we should mask/ignore when detecting phase changes (by
filtering VCFs accordingly - so this would involve updating vcfprep down the line)

found it - looks like it's in

```
/research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/misc_analyses/ \
chlamy_v6_project/QC17_CC4532/SNP_filtering
```

there's also a `QC17_v6` VCF but I recall Rory cautioning me not to use this one, given
that it has some manually inserted 'null' calls I likely won't want to deal with 

going to symlink this into `data/references/` as `QC_v6_all.vcf.gz` as well as the unfiltered
version in the folder above, since that may be useful in finding 'crap' regions

```
lrwxrwxrwx. 1 domain users  169 Nov  2 15:49 QC_v6_all.vcf.gz -> /research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/misc_analyses/chlamy_v6_project/QC17_CC4532/SNP_filtering/QC17_CC4532.all_sites.all_isolates.vcf.gz
lrwxrwxrwx. 1 domain users  173 Nov  2 15:49 QC_v6_all.vcf.gz.tbi -> /research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/misc_analyses/chlamy_v6_project/QC17_CC4532/SNP_filtering/QC17_CC4532.all_sites.all_isolates.vcf.gz.tbi
lrwxrwxrwx. 1 domain users  146 Nov  2 15:49 QC_v6_unfiltered.vcf.gz -> /research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/misc_analyses/chlamy_v6_project/QC17_CC4532/QC17_CC4532.GenotypeGVCFs.vcf.gz
lrwxrwxrwx. 1 domain users  150 Nov  2 15:49 QC_v6_unfiltered.vcf.gz.tbi -> /research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/misc_analyses/chlamy_v6_project/QC17_CC4532/QC17_CC4532.GenotypeGVCFs.vcf.gz.tbi
```

## 11/1/2022

today - regenerating all VCFs that involve 1691 since that BAM was remade

snakemake is being a total imbecile and refusing to detect that the CC1691 bam
was regenerated - not sure what's going on here - so I'm going to have to run the command
manually I suppose

```bash
time gatk HaplotypeCaller \
-R data/references/CC4532.w_organelles_MTplus.fa \
-I data/alignments/parental_bam/CC1691.bam \
-ERC GVCF -ploidy 2 \
--heterozygosity 0.02 --indel-heterozygosity 0.002 \
--native-pair-hmm-threads 2 -O data/genotyping/gvcf_sample/CC1691.g.vcf.gz
```

## 13/1/2022

alright - snakemake is still being uncooperative and not detecting that something changed -
I'm just going to remove all the final `vcf_filtered` files that involve 1691 in some way
into `temp_1691/` (in the project root) and get snakemake to go again

for some reason it now feels the need to redo all the gvcf genotyping again, which is
just great, but I might as well just let it do its thing since that step is quick enough

also updated vcfprep to _just output SNPs_ in VCFs based on last Rob meeting, since we
can revisit re: indels down the line - though this does mean I have to run all instances
of readcomb (parental and regular) once more since the marker sets have changed

other than this, I need to (use that combined v6 VCF to?) do allele depth scans and 
look for regions that are potentially paralogous - alternatively, I should add a purity
filter to `vcfprep` (which would be handy to have regardless) 

found another bug too... here's a to-do list:

- fix parental allele difference check - currently, if `gt_bases` is `['G|G', 'G/G']` that
gets counted as an informative allele
- add purity filter

## 14/1/2022

all done and implemented in 0.1.6 - here we go again

```bash
# redoes readcomb_vcfprep and vcf_tabix rules
time snakemake -pr -s analysis/genotyping/variant_calling.smk --cores 8
```



