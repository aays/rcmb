
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

## 18/10/2022

today - getting denominators for annotations 

for this, I need to iterate through the relevant GFF features
(CDS, 5' UTR, 3' UTR) and also calculate intergenic
sequence (between genes, between CDS) 

nts: I know genes sometimes overlap, but do CDS ever overlap? 
let's just assume they do - if so, there just won't be an intron there

for each feature, need to look in the `read_counts` (`data/callability/tracts`)
files to then just get the callable intervals within that feature, and get the
denominators for that

giving this a go:

```bash
time python analysis/rate/annotation_effective_sequence.py \
--bam data/alignments/bam_prepped/sorted/2343x1691.sorted.bam \
--gff data/references/CC4532.v1_1.genes.primaryTrs.db \
--snp_read_counts data/callability/tracts/read_counts_tabix/2343x1691.read_counts.tsv.gz \
--out ant_test.tsv
```

this is genuinely going to run forever - let's try this again with a simplified approach
that just throws out any features that are entirely SNP deserts - getting more granular
is going to take ages

might be worth expanding this to 'first exon', 'second exon', etc

## 19/10/2022

alright - so this tidily wrapped up in just over an hour, and looks great

I'm going to get this running over snakemake for now, and then
we can calculate the annotation rates for 2343x1691 in a terminal
just as a sanity check

```bash
time snakemake -pr -s analysis/rate/rate.smk --cores 10 # should be done in 3 hours - took 3.5
```

next up - need to get annotations for each of the COs to then
calculate the number of COs per annotation

I think I can do this in the terminal - something like:

```python
import csv
import gffutils

db = gffutils.FeatureDB('data/references/CC4532.v1_1.genes.primaryTrs.db')

with open('data/rate/crossovers_all_annotations.tsv', 'w') as f_out:
    fieldnames = [
        'cross', 'chrom', 'midpoint', 'read_name',
        'CDS', 'five_prime_UTR', 'three_prime_UTR', 'intronic', 'intergenic']
    writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    with open('data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        for line in tqdm(reader):
            midpoint = line['midpoint']
            region = midpoint + '-' + midpoint.split(':')[-1]
            
            out_dict = {
                'cross': line['cross'],
                'chrom': line['chrom'],
                'midpoint': midpoint,
                'read_name': line['read_name'],
                'CDS': 0, 'five_prime_UTR': 0, 
                'three_prime_UTR': 0, 'intronic': 0, 'intergenic': 0
                }

            features_overlap = sorted(list(set([feature.featuretype 
                for feature in db.region(region)])))
            if not features_overlap:
                out_dict['intergenic'] = 1
            elif (
                'gene' in features_overlap and 
                'CDS' not in features_overlap and
                'three_prime_UTR' not in features_overlap and
                'five_prime_UTR' not in features_overlap
            ):
                out_dict['intronic'] = 1
            for feature in features_overlap:
                if feature in out_dict:
                    out_dict[feature] += 1
            writer.writerow(out_dict)
```
