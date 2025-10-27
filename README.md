# Constellation illuminates rare disease genetics

**This is the github repository to supplement our preprint : https://www.medrxiv.org/content/10.1101/2025.10.15.25337675v1**

The repository contains scripts and documentation for assessing the performance of the BCM-Constellation variant calling pipeline. The assessment includes benchmarking based on HG002 samples and a annotation/prioritization pipeline for solved and unsolved rare disease trios.

There are two main components in this repository: HG002 Benchmarking pipeline and Annotation pipeline.

## Data Availability
All VCF files used in HG002 benchmark were uploaded to: https://zenodo.org/records/17297907

## Annotation pipeline

**annotsv_individual.sh**

Description:
Runs AnnotSV
 on individual structural variant (SV) VCF files. Each SLURM array job processes one VCF file listed in vcflist.txt. It annotates SVs against the specified annotation directory and optionally produces both .tsv and .vcf outputs.

```
Inputs:

LIST_FILE — text file listing input VCF paths (default: vcflist.txt)

ANNOTSV_DIR — AnnotSV installation directory

ANN_DIR — AnnotSV annotation database directory

GENES_FILE (optional) — candidate gene list file

EMIT_VCF (flag) — whether to output annotated VCF (1) or not (0)

Outputs:

<OUT_DIR>/<SAMPLE>/
├─ <sample>.annotsv.tsv — annotated SV table
├─ <sample>.annotsv.vcf.gz (optional) — annotated VCF
```


**annotsv_trio.sh**

Description:
Performs trio-based SV analysis using AnnotSV. Each SLURM array job processes one trio (proband–mother–father) defined in trios.tsv. The script merges the three VCFs, identifies proband-specific and inherited variants, and annotates both sets with AnnotSV.

```
Inputs:

TRIO_TSV — TSV file with 3 columns: <proband_vcf> <mother_vcf> <father_vcf>

ANNOTSV_DIR — AnnotSV installation directory

ANN_DIR — AnnotSV annotation directory

GENOME_BUILD — reference genome (default: GRCh38)

Outputs:

<WORK_DIR>/<TRIO_TAG>/
├─ trio.merged.vcf.gz — merged trio VCF
├─ proband_specific_strict.vcf.gz — proband-only variants
├─ inherited_any_parent.vcf.gz — inherited variants from either parent
├─ <TRIO_TAG>.proband_specific_strict.annotsv.tsv — AnnotSV annotation for proband-specific variants
├─ <TRIO_TAG>.inherited_any_parent.annotsv.tsv — AnnotSV annotation for inherited variants

```

**annovar_annotation_array.sh**

Description:
Annotates trio SNV/indel data using ANNOVAR. It normalizes and merges VCFs from each family trio, computes Mendelian error (MER) statistics via bcftools +mendelian2, and annotates resulting variant sets (good, mendelian error, annotated) with ANNOVAR.

```
Inputs:

TRIO_TSV — TSV with 4 columns: <proband_vcf> <mother_vcf> <father_vcf> <sex>
(Sex: 1 for male, 2 for female)

REF_FASTA — reference genome FASTA (hg38)

ANNOVAR_DIR — path to ANNOVAR installation

HUMANDB — path to ANNOVAR humandb directory

(Optional) REGION_BED — BED file to restrict region of analysis

Outputs:

annovar_out/<PROBAND_ID>/
├─ trio.merged.vcf.gz / trio.renamed.vcf.gz — merged & normalized trio VCF
├─ mendel.stats.txt — Mendelian statistics summary
├─ trio.merr.annotated.vcf.gz — Mendelian error variants
├─ trio.good.annotated.vcf.gz — Mendelian-consistent variants
├─ trio.error.vcf.gz — Mendelian error-only variants
├─ trio_merr_annot.hg38_multianno.txt — ANNOVAR annotations for MER variants
├─ trio_good_annot.hg38_multianno.txt — ANNOVAR annotations for consistent variants
├─ trio_error_annot.hg38_multianno.txt — ANNOVAR annotations for error variants

```

**bcftools_mer_slurm.sh**

Description:
Automates the computation of Mendelian error rates for trios using bcftools +mendelian2. It merges parental and proband VCFs, indexes them, and runs the plugin to extract summary statistics or custom MER TSVs via an external converter script.

```
Inputs:

TSV — TSV with 4 columns: <proband_vcf> <mother_vcf> <father_vcf> <sex>

MER_SCRIPT — path to a script that parses mendelian2 summaries into MER tables

Outputs:

<OUTDIR>/<PROBAND_ID>.merged.vcf.gz — merged VCF

<OUTDIR>/<PROBAND_ID>.merged.vcf.gz.tbi — tabix index

<OUTDIR>/<PROBAND_ID>_MER.tsv — Mendelian error summary table (via MER_SCRIPT)

```

**extract_af0.001_rm_synomousSNV_novel_in_proband.py**

Description:
Filters annotated ANNOVAR trio TSV files to extract rare (AF < 0.001), non-synonymous, proband-specific (de novo), and Mendelian error SNVs. It normalizes genotypes, detects inheritance patterns, and outputs three tables summarizing variant categories.

```
Inputs:

input — ANNOVAR trio annotation TSV (*.multianno.txt)

Parameters:

--af-col — column for allele frequency (default: gnomad41_genome_AF)

--af-threshold — AF cutoff (default: 0.001)

--format-col, --proband-col, --mother-col, --father-col — optional manual overrides

Outputs:

<prefix>_altcarriers.tsv — all variants where proband carries an ALT allele

<prefix>_denovo.tsv — strict de novo variants (absent in both parents)

<prefix>_mendelian_errors.tsv — variants violating Mendelian inheritance

```

**extract_het_compound_snv.py**

Description:
Identifies compound heterozygous SNVs from ANNOVAR trio annotations. The script detects cases where a proband has two heterozygous variants in the same gene, inherited from different parents (one maternal, one paternal), while filtering out synonymous and common variants.

```
Inputs:

input — ANNOVAR trio TSV (*.multianno.txt)

Parameters:

--af-col — allele frequency column (default: gnomad41_genome_AF)

--af-threshold — maximum AF threshold (default: 0.001)

--format-col, --proband-col, --mother-col, --father-col — optional overrides

Outputs:

<prefix>_comphet_lenient.tsv — variants forming potential compound-heterozygous pairs
(includes both maternal- and paternal-inherited candidates)
```

## HG002 Benchmarking pipeline

**bwa_pe.sh** – Aligns paired-end FASTQ reads to a reference genome using BWA, then sorts and indexes the resulting BAM file.

**clair3.sh** – Runs Clair3 variant calling on an input BAM file with specified models.

Models: <r1041_e82_400bps_sup_v520 for ONT> <hifi_revio for PB>

**deepvariant_wgs.sh** – Run DeepVariant to perform whole-genome (WGS) variant calling in Illumina standard short-read DNA sequncing datasets.

**manta.sh** – Detects structural variants (SVs) from Illumina standard short-read DNA sequncing datasets.

**rtgtools_vcfeval.sh** – Evaluates small variant calling accuracy.

**sniffles2_benchmark.sh** – Runs Sniffles2 for structural variant detection in long-read datasets(PacBio and ONT)

**truvariv5.3.sh** – Benchmarks SV callsets using Truvari.

**whatshap_phase.sh** – Performs haplotype phasing of variants in callsets except Constellation.

**whatshap_stats.sh**  Generates phasing statistics and block summaries from a phased VCF, this script was performed on all callsets including the Constellation callset.


## Publication
Check out the preprint [here]().
