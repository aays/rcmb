
## 7/10/2022

today - motif searches in COs using meme 5.4.1, which has been installed
in `~/apps` and can be found in `~/meme/bin`

looks like MEME requires a fasta containing sequences of interest - let's do that
for the three populations individually, and then all combined

based on my brief MEME notebook from late 2017 - I can use `bedtools getfasta` to generate this
fasta file:

```bash
# first need to get bed files from the CO files

mkdir -p data/motifs/co_bed/

for fname in data/phase_changes/crossovers_lasso/*; do
    base=$(basename ${fname})
    grep -v '3071x2931' ${fname} | grep -v '3071x3062' | cut -f 5,7,8 | tail -n +2 > data/motifs/co_bed/${base}
done

# sort - within co_bed dir
for fname in *.tsv; do
    base=$(basename ${fname} .tsv)
    sort -k1,1 -k2,2n ${fname} > ${base}.sorted.tsv;
done

# bedtools getfasta - back out to project dir
mkdir -p data/motifs/fasta

for fname in data/motifs/co_bed/*sorted.tsv; do
    base=$(basename ${fname} .sorted.tsv)
    bedtools getfasta \
        -fi data/references/CC4532.w_organelles_MTplus.fasta \
        -bed ${fname} \
        -fo data/motifs/fasta/${base}.fasta
done
```

and now let's get MEMEing - 

```bash
mkdir -p data/motifs/meme_out

# test on smallest file first
time ./bin/meme data/motifs/fasta/crossovers_NA1_corrected.fasta \
-o data/motifs/meme_out/ \
-text -dna -V -p 8 \
-nmotifs 10 \ # find up to 10 motifs
-minsites 5 -maxsites 20
```

this is just a preliminary run - realistically I should be using the actual
sequence of the parents/reads for this, with each the sequence of each rcmb event
considered + extended about ~500 bp in either direction (Comeron 2012)

the 'correct' way to do this would otherwise be to, for each CO:

- get parents, midpoint, and read in 'detection' list with `eval`
- get 250 bp on either side of midpoint from either parental fasta
    - eg if detection is 1,2 and midpoint is 2000, get 1750-2000 from 1 and 2000-2250 from 2
    - could use `SeqIO.index` for quick lazy loading of records this way
    - I do think I already have chrom-separated fasta files from the ldhelmet analysis though which works
- append to fasta file of CO sequences

also need to similarly sample n (number of COs) randomly sampled 500 bp chunks of the genome
and provide these to MEME with the `-neg` arg (file containing control sequences) 

trying out a script that does this:

```bash
time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_NA1_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--out data/motifs/fasta/NA1_spliced.fasta # takes 1 min per 1k COs
```

looks good - and now to run MEME on this - this does require
a custom alphabet since othe fastas have '\*' chars for some reason,
so I've copied the DNA alphabet definition and added '\*' at the bottom
as equivalent to N

```bash
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-o data/motifs/meme_out/ \
-text -V -p 8 -alph meme_alphabet.txt \
-nmotifs 10 \ # find up to 10 motifs
-minsites 5 -maxsites 20 > data/motifs/meme_out/NA1.meme
```

next up install ceqlogo and convert meme output - try generating the html
output since maybe that'll automatically generate the seq logo figure

## 8/10/2022

let's try meme html output

```bash
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-o data/motifs/meme_out_NA1/ \ # have to give new dir
-V -p 8 -alph meme_alphabet.txt \
-nmotifs 1  -minsites 5 -maxsites 20
```

this actually worked and generated the files + figures I needed! 

let's rerun with 8 motifs - once I'm convinced this looks good,
going to generate a set of randomly drawn control sequences as well

```bash
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-oc data/motifs/meme_out_NA1/ \ # will overwrite prev dir
-V -p 8 -alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20 # took 18 min
```

looks good, and this seems to be returning a single significant motif - the final
thing left to do is redo this with control sequences, and then run it on
all three populations 

updated `generate_fastas.py` so that it'll randomly draw midpoints
if needed - 

```bash
time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_NA1_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--random \
--out data/motifs/fasta/NA1_random.fasta
```

and let's try MEME one last time for this before doing this across pops:

```bash
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-neg data/motifs/fasta/NA1_random.fasta \ # randomly drawn sequences
-oc data/motifs/meme_out_NA1/ \ # will overwrite prev dir
-objfun se \ # selective enrichment
-V --alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20 # took 18 min
```

going to also try using STREME after this - apparently that is specifically
discriminative, and better for larger datasets (e.g greater than 50 sequences,
and we have thousands) 

```bash
time ./bin/streme \
-p data/motifs/fasta/NA1_spliced.fasta \
-n data/motifs/fasta/NA1_random.fasta \
-o data/motifs/streme_out_NA1 \
--alph meme_alphabet.txt \
-nmotifs 8 \
-minw 5 -maxw 20
```

there's a really convincing motif here! let's regen this for the other pops (NA2 and mix)

```bash
time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_NA2_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--out data/motifs/fasta/NA2_spliced.fasta

time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_NA2_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--random \
--out data/motifs/fasta/NA2_random.fasta

time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_mix_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--out data/motifs/fasta/mix_spliced.fasta

time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_mix_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--random \
--out data/motifs/fasta/mix_random.fasta
```

## 9/10/2022

today - running STREME on these other two populations 

```bash
time ./bin/streme -p data/motifs/fasta/NA2_spliced.fasta \
-n data/motifs/fasta/NA2_random.fasta \
-oc data/motifs/streme_out_NA2 \
--alph meme_alphabet.txt \
-nmotifs 8 \
-minw 5 \
-maxw 20 # took 10 min

time ./bin/streme -p data/motifs/fasta/mix_spliced.fasta \
-n data/motifs/fasta/mix_random.fasta \
-oc data/motifs/streme_out_mix \
--alph meme_alphabet.txt \
-nmotifs 8 \
-minw 5 \
-maxw 20
```

and finally, doing this for all crossovers combined:

```bash
time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--out data/motifs/fasta/all_spliced.fasta

time python analysis/motifs/generate_fastas.py \
--fname data/phase_changes/crossovers_lasso/crossovers_all_corrected.tsv \
--window_size 250 \
--fasta_dir data/ldhelmet/fasta \
--random \
--out data/motifs/fasta/all_random.fasta

# and then streme
time ./bin/streme -p data/motifs/fasta/all_spliced.fasta \
-n data/motifs/fasta/all_random.fasta \
-oc data/motifs/streme_out_mix \
--alph meme_alphabet.txt \
-nmotifs 8 \
-minw 5 \
-maxw 20 # took 52 min! 
```

finally - I'm going to run MEME on the full file in case,
since with the full CO set there are so many COs that
my randomly drawn sequences might actually overlap with them -
so here goes:

```bash
# running without randomly drawn sequences - since that's the STREME results
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-oc data/motifs/meme_out_NA1/ \ # will overwrite prev dir
-V -alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/NA2_spliced.fasta \
-oc data/motifs/meme_out_NA2/ \
-V -alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/mix_spliced.fasta \
-oc data/motifs/meme_out_mix/ \
-V -alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/all_spliced.fasta \
-oc data/motifs/meme_out_all/ \
-V -alph meme_alphabet.txt \
-nmotifs 8 -minsites 5 -maxsites 20
```

gonna try the central distance algorithm, which scores motifs
based on distance from the centre

```bash
# running without randomly drawn sequences - since that's the STREME results
time ./bin/meme data/motifs/fasta/NA1_spliced.fasta \
-oc data/motifs/meme_out_NA1_cd/ \ # will overwrite prev dir
-V -alph meme_alphabet.txt \
-objfun cd -nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/NA2_spliced.fasta \
-oc data/motifs/meme_out_NA2_cd/ \
-V -alph meme_alphabet.txt \
-objfun cd -nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/mix_spliced.fasta \
-oc data/motifs/meme_out_mix_cd/ \
-V -alph meme_alphabet.txt \
-objfun cd -nmotifs 8 -minsites 5 -maxsites 20 # 18 min

time ./bin/meme data/motifs/fasta/all_spliced.fasta \
-oc data/motifs/meme_out_all_cd/ \
-V -alph meme_alphabet.txt \
-objfun cd -nmotifs 8 -minsites 5 -maxsites 20
```

