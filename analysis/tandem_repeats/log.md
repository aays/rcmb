
## 4/9/2022

today - running ULTRA to get tandem repeats across all chromosomes 
using the LDhelmet per-chr fastas

```bash
mkdir -p data/tandem_repeats
```

I figured out how to run this with singularity in the `phase_changes` log - now
to parametrize it properly and scale it up here

after that, will need to convert the output jsons into something that can
be tabixed so that `summarise_cross` can look that up for a given pair

it does occur to me that the final files should be sample-specific and contain
all chrs for each sample - something to keep in mind

here's the ULTRA command that worked (though absolute paths being needed is annoying)

```bash
time singularity exec --bind \
/research/projects/chlamydomonas/genomewide_recombination/rcmb/data/ldhelmet/fasta_chrom/:\
/research/projects/chlamydomonas/genomewide_recombination/rcmb/analysis/tandem_repeats/mnt \
analysis/tandem_repeats/ultra.sif \
ultra -f test.json -at 0.36 \
-s 9 \ # score threshold
data/ldhelmet/fasta_chrom/NA2.chromosome_01.fa
```

first pass at making the json files for each chr -

```bash
time snakemake -pr -s analysis/tandem_repeats/tandem_repeats_calc.smk --cores 8
```

in the meantime - going to make a script that converts these json files to
tsv files and then add it in to the workflow

took 45 min - and now for json to tsv

```bash
time snakemake -pr -s analysis/tandem_repeats/tandem_repeats_calc.smk --cores 4
```

looks good - finally, bgzip and tabix as bed files:

```
time snakemake -pr -s analysis/tandem_repeats/tandem_repeats_calc.smk
```

