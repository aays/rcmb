# to be run from project root
if (test "$(basename ${PWD})" != "rcmb"); then
    echo "[rcmb] script needs to be run from project root"
    echo "[rcmb] aborting"
    exit 0
fi

mkdir -pv data/references
mkdir -pv data/alignments/bam
mkdir -pv data/alignments/fastq
mkdir -pv data/alignments/fastq/fastqc
mkdir -pv data/alignments/fastq/fastqc/raw
mkdir -pv data/alignments/fastq_trim
mkdir -pv data/alignments/fastq_trim/fastqc
mkdir -pv data/alignments/fastq_trim/fastqc/raw
