# to be run from project root
if (test "$(basename ${PWD})" != "rcmb"); then
    echo "[rcmb] script needs to be run from project root"
    echo "[rcmb] aborting"
    exit 0
fi

mkdir -pv data/references
mkdir -pv data/alignments/bam

for dirname in fastq fastq_trim parental_fastq parental_fastq_trim; do
    mkdir -pv data/alignments/${dirname}
    mkdir -pv data/alignments/${dirname}/fastqc
    mkdir -pv data/alignments/${dirname}/fastqc/raw
done
