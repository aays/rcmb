# `genotyping`

Genotyping workflow for parental samples.

### `variant_calling.smk`

Snakemake workflow for parental pair VCF generation. This script
follows the [GATK germline short variant discovery workflow][gatk-workflow]
and calls individual GVCFs before combining as needed and performing
joint genotyping with HaplotypeCaller.

CombineGVCFs was used instead of GenomicsDBImport since no combined GVCF had
more than two samples. The workflow also does not include Variant Quality Score
Recalibration due to the lack of a high quality reference dataset in _C.
reinhardtii_.

`readcomb-vcfprep` was used for downstream hard filtering of VCFs for high quality
variants. 

[gatk-workflow]: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
