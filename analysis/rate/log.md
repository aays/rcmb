
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


