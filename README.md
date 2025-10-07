# BCM-Constellation Assessment

## HG002 Benchmarking pipeline



## Annotation pipeline

**bwa_pe.sh** – Aligns paired-end FASTQ reads to a reference genome using BWA, then sorts and indexes the resulting BAM file.

**clair3.sh** – Runs Clair3 variant calling on an input BAM file with specified models.

Models: <r1041_e82_400bps_sup_v520 for ONT> <hifi_revio for PB>

**deepvariant_wgs.sh** – Run DeepVariant to perform whole-genome (WGS) variant calling in Illumina standard short-read DNA sequncing datasets.

**manta.sh** – Detects structural variants (SVs) from Illumina standard short-read DNA sequncing datasets.

**rtgtools_vcfeval.sh** – Evaluates small variant calling accuracy.

**sniffles2_benchmark.sh** – Runs Sniffles2 for structural variant detection in long-read datasets(PacBio and ONT)

**truvariv5.3.sh** – Benchmarks SV callsets.

**whatshap_phase.sh** – Performs haplotype phasing of variants in callsets except Constellation.

**whatshap_stats.sh**  Generates phasing statistics and block summaries from a phased VCF, this script was performed on all callsets including the Constellation callset.