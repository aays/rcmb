# `alignment`

Scripts and workflows to generate BAM files for both recombinant and parental samples.

Parental samples are handled separately since these FASTQs were obtained from multiple
previous studies (see `log.md`) whereas recombinant samples are all derived from the same
sequencing run.

### `dir_setup.sh`

Generates directory tree for both recombinant and parental sample alignment.

```bash
bash analysis/alignment/dir_setup.sh
```

### `read_qc_trim.smk` and `parental_fastq_qc.smk`

Snakemake workflows to QC reads, trim with trimmomatic, and then QC post-trimming.

These workflows are meant to be run one rule at a time, with trimming parameters
adjusted as needed based on QC results.

```bash
# e.g.
snakemake -pr -s read_qc_trim.smk fastqc # run fastqc on recombinant samples
```

### `alignment.smk` and `parental_alignment.smk`

Full alignment workflows (FASTQ -> BAM) for both recombinant and parental samples. 

```bash
snakemake -pr -s alignment.smk --cores 4
snakemake -pr -s parental_alignment.smk --cores 4
```
