
## 3/10/2022

main two analyses here -

1. how distant are COs from their nearest TSS? 
2. how distant are COs from H3K4me3 peaks?

TSS data - need to use the GFF for - and peaks I've obtained
from Rory (bed file named `H3K4me3_peaks` in narrowPeak format)

for TSS - will need to create some kind of bed file of TSS's from
the gff and do quick lookups for each CO to get its distance
to its nearest TSS

for chromatin marks - already have the bed file, just need
to get distance from midpoint

## 4/10/2022

let's start with TSS, since it actually might be the more annoying of the two

going to fetch all the mRNA records in the gff and create a bed file out of that

```python
import csv
import gffutils
from tqdm import tqdm

db = gffutils.FeatureDB('data/references/CC4532.v1_1.gene_exons.db', keep_order=True)

with open('data/tss_marks/tss_sites.tsv', 'w') as f:
    fieldnames = ['#chromosome', 'tss_pos', 'strand', 'start', 'end']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for i, feature in tqdm(enumerate(db.all_features(featuretype='mRNA'))):
        out_dict = {
            '#chromosome': feature.chrom, 'start': feature.start', 
            'end': feature.end, 'strand': feature.strand}
        if feature.strand == '+':
            out_dict['tss_pos'] = feature.start
        elif feature_strand == '-':
            out_dict['tss_pos'] = feature.end
        writer.writerow(out_dict)
```

next up - for each CO, need to calculate distance to nearest TSS 

let's zip and tabix this after sorting - 

```bash
sort -k1,1 -k2,2n tss_sites.tsv > tss_sites_sorted.tsv

bgzip tss_sites_sorted.tsv
tabix -p vcf tss_sites_sorted.tsv.gz
```

and now let's go through the cos and get proximity for each - I could
make this a reusable script but I don't have a ton of time so _shrug_ 

```python
import math
import csv
import pysam
from tqdm import tqdm

tss = pysam.TabixFile('data/tss_marks/tss_sites_sorted.tsv.gz')

def get_tss_dist(chrom, co_pos):
    tss_chrom = tss.fetch(chrom)
    min_pos = min([abs(int(line.split('\t')[1]) - co_pos) for line in tss_chrom])
    # re-instantiate generator
    tss_chrom = tss.fetch(chrom)
    nearest_tss = [line for line in tss_chrom if abs(int(line.split('\t')[1]) - co_pos) == min_pos][0]
    _, tss_pos, strand, _, _ = nearest_tss.split('\t')
    tss_pos = int(tss_pos)
    if strand == '+':
        dist = co_pos - tss_pos
    elif strand == '-':
        dist = -1 * (co_pos - tss_pos)
    return dist, tss_pos
        

with open('data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv', 'r') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    with open('data/tss_marks/co_tss_distances.tsv','w') as f_out:
        fieldnames = [
            'cross', 'type', 'chromosome', 'midpoint', 
            'start', 'end', 'read_name', 'tss_nearest', 'tss_dist']
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader() 
        for line in tqdm(reader):
            chrom = line['chromosome']
            midpoint = int(line['midpoint'].split(':')[1])
            dist, tss_pos = get_tss_dist(chrom, midpoint)

            out_dict = {
                'cross': line['cross'],
                'type': line['type'],
                'chromosome': line['chromosome'],
                'midpoint': line['midpoint'],
                'start': line['start'],
                'end': line['end'],
                'read_name': line['read_name'],
                'tss_nearest': tss_pos,
                'tss_dist': dist
                }
            writer.writerow(out_dict)

```

and now doing the same for the H3K4 peaks - will need to get distance from centre of peak,
however

first need to bgzip and tabix the narrowPeaks file to make it tabix-usable

```bash
# in data/tss_marks
bgzip H3K4me3_peaks.narrowPeak
tabix -p bed H3K4me3_peaks.narrowPeak.gz
```

looks good - and now to generate distances:

```python
import math
import csv
import pysam
from tqdm import tqdm

peaks = pysam.TabixFile('data/tss_marks/H3K4me3_peaks.narrowPeak.gz')

def get_peak_dist(chrom, co_pos):
    peak_chrom = peaks.fetch(chrom)
    peak_midpoints_all = []
    for line in peak_chrom:
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        peak_midpoint = int(start + (end - start))
        peak_midpoints_all.append(peak_midpoint)
    min_dist = min([abs(peak_midpoint - co_pos) for peak_midpoint in peak_midpoints_all])
    peak_chrom = peaks.fetch(chrom)
    for line in peak_chrom:
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        peak_midpoint = int(start + (end - start))
        if abs(peak_midpoint - co_pos) == min_dist:
            nearest_peak = line
            break
    chrom, start, end = nearest_peak.split('\t')[:3]
    return start, end

with open('data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv', 'r') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    with open('data/tss_marks/co_h3k4me3_distances.tsv','w') as f_out:
        fieldnames = [
            'cross', 'type', 'chromosome', 'midpoint', 'start', 'end', 
            'read_name', 'peak_nearest_start', 'peak_nearest_end',
            'peak_nearest_midpoint', 'peak_dist']
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader() 
        for line in tqdm(reader):
            chrom = line['chromosome']
            if chrom in ['cpDNA', 'mtDNA']:
                continue
            midpoint = int(line['midpoint'].split(':')[1])
            peak_start, peak_end = get_peak_dist(chrom, midpoint)

            out_dict = {
                'cross': line['cross'],
                'type': line['type'],
                'chromosome': line['chromosome'],
                'midpoint': line['midpoint'],
                'start': line['start'],
                'end': line['end'],
                'read_name': line['read_name'],
                'peak_nearest_start': peak_start,
                'peak_nearest_end': peak_end,
                'peak_nearest_midpoint': int(peak_start) + int((int(peak_end) - int(peak_start)) / 2)
                }
            out_dict['peak_dist'] = midpoint - out_dict['peak_nearest_midpoint']
            writer.writerow(out_dict)
```

## 5/10/2022

alright - I need to take more into account than just CO counts, which I think means
needing a script

this script could handle both peaks and TSS since there's not much difference in getting
TSS + midpoints of peaks - just need different functions to process each

other things needed would be the effective sequence denominator (for rates),
the count of SNPs, SNP density in window sizes, and parental coverage

I only have parental coverage in 2 kb overlapping windows, but that should be
fine for an approx understanding of how callable the region is (this is more of an
on/off filter I think) 

also need to report the CO score from the model! 

going to generate samtools sorted versions of the bamprep BAMs for quick denominator calc -

```bash
mkdir -p data/alignments/bam_prepped/sorted

# in bam_prepped
time for fname in *bam; do
    echo ${fname}
    samtools sort -@8 ${fname} > sorted/${fname}
done # took 1h 20 min

# and then samtools index each of them in bam_prepped/sorted
for fname in *bam; do
    echo ${fname}
    samtools index ${fname}
done
```

alright, this took all morning but let's give it a go:

```bash
time python analysis/tss_marks/get_distance.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--tss data/tss_marks/tss_sites_sorted.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/2343x1691.sorted.bam \
--vcf data/genotyping/vcf_filtered/2343x1691.vcf.gz \
--window_size 100 \
--out data/tss_marks/tss_marks_all.tsv
```

## 6/10/2022

generating a neutral expectation for the COs - need a script that will
take counts of COs per chr per sample and then output a file similar
to my all COs file except with evenly spaced positions 

giving this a go:

```bash
time python analysis/tss_marks/neutral_cos.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--out data/tss_marks/neutral_cos.tsv
```

looks good - and now to get the TSS/marks distances:

WAIT - there's a bug - I've been using the 2343x1691 bam and VCF for _all_
of the windows! let's fix up the script so that each cross gets its own
proper denominator and SNP density

```bash
time python analysis/tss_marks/get_distance.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--tss data/tss_marks/tss_sites_sorted.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_all.tsv
```

and then also running this on the neutral COs anyways:

```bash
time python analysis/tss_marks/get_distance.py \
--fname data/tss_marks/neutral_cos.tsv \
--tss data/tss_marks/tss_sites_sorted.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_neutral.tsv
```

## 10/10/2022

I still did this wrong - these `eff_bp` denominators were wrong again - 
let's rerun this using the `truncate` option, which should calculate
the bases in a given window correctly 

```bash
# also updated the script to ignore the contaminated crosses
# regular COs
time python analysis/tss_marks/get_distance.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--tss data/tss_marks/tss_sites_sorted.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_all.tsv

# randomly drawn COs
time python analysis/tss_marks/get_distance.py \
--fname data/tss_marks/neutral_cos.tsv \
--tss data/tss_marks/tss_sites_sorted.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_neutral.tsv
```

## 18/10/2022

I think I've been using the wrong GFF! just used the `primaryTrs` GFF and there
are far fewer (17k instead of 32k) TSSes - maybe this is why things have been going wrong

going to move the alt one into our data folder and do this again:

```bash
# generated TSS file using code at the top of this log, just with other GFF
# in data/tss_marks/

sort -k1,1 -k2,2n tss_sites_alt.tsv > tss_sites_sorted_alt.tsv
bgzip tss_sites_sorted_alt.tsv
tabix -p vcf tss_sites_sorted_alt.tsv.gz

time python analysis/tss_marks/get_distance.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--tss data/tss_marks/tss_sites_sorted_alt.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_all_alt.tsv

# randomly drawn COs
time python analysis/tss_marks/get_distance.py \
--fname data/tss_marks/neutral_cos.tsv \
--tss data/tss_marks/tss_sites_sorted_alt.tsv.gz \
--peaks data/tss_marks/H3K4me3_peaks.narrowPeak.gz \
--bam data/alignments/bam_prepped/sorted/ \
--vcf data/genotyping/vcf_filtered/ \
--window_size 100 \
--out data/tss_marks/tss_marks_neutral_alt.tsv
```
