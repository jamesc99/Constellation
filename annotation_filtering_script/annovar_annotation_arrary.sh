#!/bin/bash
#SBATCH --job-name=trio_annovar
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

TRIO_TSV="${1:?Usage: sbatch --array=1-N annovar_annotation_array.sh <trios.tsv>}"

REF_FASTA="${REF_FASTA:?Set REF_FASTA to hg38 FASTA path}"
ANNOVAR_DIR="${ANNOVAR_DIR:?Set ANNOVAR_DIR to the annovar install dir}"
HUMANDB="${HUMANDB:?Set HUMANDB to the annovar humandb dir}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
TABIX_BIN="${TABIX_BIN:-tabix}"
THREADS="${THREADS:-8}"
ANNOVAR_THREADS="${ANNOVAR_THREADS:-8}"

PROTOCOLS="${PROTOCOLS:-refGene,refGeneWithVer,avsnp151,clinvar_20240611,gnomad41_genome,gnomad211_exome,dbnsfp47a,dbscsnv11,cytoBand}"
OPERATIONS="${OPERATIONS:-g,g,f,f,f,f,f,f,r}"

REGION_BED="${REGION_BED:-}"  # optional, leave empty to use whole genome

read -r PROBAND_VCF MOTHER_VCF FATHER_VCF PROBAND_SEX < <(
  awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1"\t"$2"\t"$3"\t"$4}' "${TRIO_TSV}"
)

PROBAND_ID="$(basename "${PROBAND_VCF}" | sed 's/_corrected\.hard-filtered\.vcf\.gz$//')"
MOTHER_ID="$(basename  "${MOTHER_VCF}"  | sed 's/_corrected\.hard-filtered\.vcf\.gz$//')"
FATHER_ID="$(basename  "${FATHER_VCF}"  | sed 's/_corrected\.hard-filtered\.vcf\.gz$//')"

if [[ "${PROBAND_SEX}" == "1" ]]; then
  PFM="1X:${PROBAND_ID},${FATHER_ID},${MOTHER_ID}"
elif [[ "${PROBAND_SEX}" == "2" ]]; then
  PFM="2X:${PROBAND_ID},${FATHER_ID},${MOTHER_ID}"
else
  echo "[ERROR] proband_sex must be 1 (male) or 2 (female)"; exit 1
fi

[[ -s "${PROBAND_VCF}" && -s "${MOTHER_VCF}" && -s "${FATHER_VCF}" ]] || { echo "[ERROR] missing input VCF(s)"; exit 1; }
[[ -s "${REF_FASTA}" && -s "${REF_FASTA}.fai" ]] || { echo "[ERROR] missing REF_FASTA or .fai"; exit 1; }
[[ -d "${HUMANDB}" ]] || { echo "[ERROR] humandb not found: ${HUMANDB}"; exit 1; }

BASE_OUTDIR="annovar_out/${PROBAND_ID}"
mkdir -p "${BASE_OUTDIR}"

declare -A MAP=( ["proband"]="${PROBAND_VCF}" ["mother"]="${MOTHER_VCF}" ["father"]="${FATHER_VCF}" )
for role in proband mother father; do
  in_vcf="${MAP[$role]}"
  out_vcf="${BASE_OUTDIR}/${role}.norm.vcf.gz"
  "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -f "${REF_FASTA}" -m -any -Oz -o "${out_vcf}" "${in_vcf}"
  "${TABIX_BIN}" -f "${out_vcf}"
done

"${BCFTOOLS_BIN}" merge --threads "${THREADS}" \
  -Oz -o "${BASE_OUTDIR}/trio.merged.vcf.gz" \
  "${BASE_OUTDIR}/proband.norm.vcf.gz" \
  "${BASE_OUTDIR}/mother.norm.vcf.gz" \
  "${BASE_OUTDIR}/father.norm.vcf.gz"
"${TABIX_BIN}" -f "${BASE_OUTDIR}/trio.merged.vcf.gz"

printf "%s\n%s\n%s\n" "${PROBAND_ID}" "${MOTHER_ID}" "${FATHER_ID}" > "${BASE_OUTDIR}/sample_order.txt"
"${BCFTOOLS_BIN}" reheader -s "${BASE_OUTDIR}/sample_order.txt" "${BASE_OUTDIR}/trio.merged.vcf.gz" \
  | "${BCFTOOLS_BIN}" view -Oz -o "${BASE_OUTDIR}/trio.renamed.vcf.gz"
"${TABIX_BIN}" -f "${BASE_OUTDIR}/trio.renamed.vcf.gz"
TRIO_VCF="${BASE_OUTDIR}/trio.renamed.vcf.gz"

MENDEL_BASE_ARGS=( +mendelian2 -p "${PFM}" )
if [[ -n "${REGION_BED}" ]]; then
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m c -R "${REGION_BED}" "${TRIO_VCF}" > "${BASE_OUTDIR}/mendel.stats.txt"
else
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m c "${TRIO_VCF}" > "${BASE_OUTDIR}/mendel.stats.txt"
fi

if [[ -n "${REGION_BED}" ]]; then
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m a -R "${REGION_BED}" -Oz -o "${BASE_OUTDIR}/trio.merr.annotated.vcf.gz" "${TRIO_VCF}"
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m ag -R "${REGION_BED}" -Oz -o "${BASE_OUTDIR}/trio.good.annotated.vcf.gz" "${TRIO_VCF}"
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m e -R "${REGION_BED}" -Oz -o "${BASE_OUTDIR}/trio.error.vcf.gz" "${TRIO_VCF}"
else
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m a  -Oz -o "${BASE_OUTDIR}/trio.merr.annotated.vcf.gz"  "${TRIO_VCF}"
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m ag -Oz -o "${BASE_OUTDIR}/trio.good.annotated.vcf.gz"  "${TRIO_VCF}"
  "${BCFTOOLS_BIN}" "${MENDEL_BASE_ARGS[@]}" -m e  -Oz -o "${BASE_OUTDIR}/trio.error.vcf.gz"         "${TRIO_VCF}"
fi
"${TABIX_BIN}" -f "${BASE_OUTDIR}/trio.merr.annotated.vcf.gz"
"${TABIX_BIN}" -f "${BASE_OUTDIR}/trio.good.annotated.vcf.gz"
"${TABIX_BIN}" -f "${BASE_OUTDIR}/trio.error.vcf.gz"

for label in merr.annotated good.annotated error; do
  IN_VCF="${BASE_OUTDIR}/trio.${label}.vcf.gz"
  OUT_PREFIX="${BASE_OUTDIR}/trio_${label}_annot"
  perl "${ANNOVAR_DIR}/table_annovar.pl" --thread "${ANNOVAR_THREADS}" \
    "${IN_VCF}" "${HUMANDB}" \
    -buildver hg38 -out "${OUT_PREFIX}" -remove -vcfinput \
    -protocol "${PROTOCOLS}" \
    -operation "${OPERATIONS}" -nastring . -polish -otherinfo
done

echo "[DONE] ${BASE_OUTDIR}"

