
## 26/9/2022

today - getting windowed denominators across all crosses using
`windowed_effective_sequence.py` (see the `phase_changes` log) 

putting together a workflow and generating 2 kb windows to start with:

```bash
time snakemake -pr -s analysis/rate/rate.smk --cores 12
# took 4 hours
```

## 27/9/2022

next up - need to get SNP counts in similar windows, simply because
SNP density is a super easy measure of divergence for the time being

this requires a much simpler script - going to just make a snakemake helper
script instead of something more complicated

giving this a go:

```bash
time snakemake -pr -s analysis/rate/rate.smk --cores 12
```

## 11/10/2022

going to redo this based on the callability (see log) updates - I now have
inter-SNP tracts and I know which ones are deserts - need to just get
denominators for windows/window chunks that are not in deserts

this will iterate through the read counts files in callability and
use that to get pileups for windows - but for a given window,
if it overlaps a SNP desert, only sequence up to the start of the SNP
desert is counted

alternatively, if the region is spanned by a SNP desert,
then the region just isn't counted at all

updating `windowed_effective_sequence.py` - I need to iterate
through the SNP tracts and then just determine windows around that

needed to generate tabixed versions in `data/callability/tracts/read_counts_tabix`
by removing the cross column and adding a '#' to the header - testing:

```bash
time python analysis/rate/windowed_effective_sequence.py \
--bam data/alignments/bam_prepped/sorted/2343x1691.sorted.bam \
--snp_read_counts data/callability/tracts/read_counts_tabix/2343x1691.read_counts.tsv.gz \
--window_size 1000 --out eff_bp_test.tsv
```
