#!/usr/bin/env bash
#SBATCH --job-name=trio_annotsv
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=medium
# Submit with: sbatch --array=1-N annotsv_trio_proband_specific_and_inherited_parent.sh <trios.tsv> [work_dir]

set -euo pipefail

TRIO_TSV="${1:?Usage: $0 <trios.tsv> [work_dir]   # TSV with 3 cols: proband_vcf<TAB>mother_vcf<TAB>father_vcf}"
WORK_DIR="${2:-trio_mode_out}"

ANNOTSV_DIR="${ANNOTSV_DIR:?Set ANNOTSV_DIR to AnnotSV install dir}"
ANNOTSV_BIN="${ANNOTSV_BIN:-${ANNOTSV_DIR}/bin/AnnotSV}"
ANN_DIR="${ANN_DIR:?Set ANN_DIR to AnnotSV annotations directory}"
GENOME_BUILD="${GENOME_BUILD:-GRCh38}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not in PATH" >&2; exit 1; }
command -v bgzip    >/dev/null 2>&1 || { echo "ERROR: bgzip not in PATH" >&2; exit 1; }
command -v tabix    >/dev/null 2>&1 || { echo "ERROR: tabix not in PATH" >&2; exit 1; }
[[ -x "${ANNOTSV_BIN}" ]] || { echo "ERROR: AnnotSV not found: ${ANNOTSV_BIN}" >&2; exit 1; }
[[ -d "${ANN_DIR}" ]] || { echo "ERROR: annotations dir not found: ${ANN_DIR}" >&2; exit 1; }
[[ -f "${TRIO_TSV}" ]] || { echo "ERROR: TRIO_TSV not found: ${TRIO_TSV}" >&2; exit 1; }

N_LINES="$(wc -l < "${TRIO_TSV}")"
[[ "${SLURM_ARRAY_TASK_ID:-1}" -ge 1 && "${SLURM_ARRAY_TASK_ID}" -le "${N_LINES}" ]] || {
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-unset} out of range (1..${N_LINES})" >&2; exit 1; }

LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TRIO_TSV}" | tr -d '\r')"
PRO_VCF="$(awk -F'\t' '{print $1}' <<< "${LINE}")"
MOM_VCF="$(awk -F'\t' '{print $2}' <<< "${LINE}")"
DAD_VCF="$(awk -F'\t' '{print $3}' <<< "${LINE}")"

for f in "${PRO_VCF}" "${MOM_VCF}" "${DAD_VCF}"; do
  [[ -f "${f}" ]] || { echo "ERROR: missing VCF: ${f}" >&2; exit 1; }
done

mkdir -p "${WORK_DIR}"

smart_index () {
  local v="$1"
  if [[ "${v}" =~ \.vcf\.gz$ ]]; then
    [[ -f "${v}.tbi" || -f "${v}.csi" ]] || tabix -f -p vcf "${v}"
  elif [[ "${v}" =~ \.vcf$ ]]; then
    bgzip -f "${v}"
    tabix -f -p vcf "${v}.gz"
    v="${v}.gz"
  else
    echo "ERROR: unsupported VCF extension: ${v}" >&2; return 1
  fi
}

smart_index "${PRO_VCF}"
smart_index "${MOM_VCF}"
smart_index "${DAD_VCF}"

PRO_ID="$(bcftools query -l "${PRO_VCF}" | head -n1)"
MOM_ID="$(bcftools query -l "${MOM_VCF}" | head -n1)"
DAD_ID="$(bcftools query -l "${DAD_VCF}" | head -n1)"
[[ -n "${PRO_ID}" && -n "${MOM_ID}" && -n "${DAD_ID}" ]] || { echo "ERROR: could not read sample IDs from VCF headers" >&2; exit 1; }

TRIO_TAG="${PRO_ID}"
TRIO_DIR="${WORK_DIR}/${TRIO_TAG}"
mkdir -p "${TRIO_DIR}"

MERGED_VCF="${TRIO_DIR}/trio.merged.vcf.gz"
bcftools merge --threads "${THREADS}" -m none -Oz -o "${MERGED_VCF}" "${PRO_VCF}" "${MOM_VCF}" "${DAD_VCF}"
tabix -f -p vcf "${MERGED_VCF}"

PS_VCF="${TRIO_DIR}/proband_specific_strict.vcf.gz"
INH_VCF="${TRIO_DIR}/inherited_any_parent.vcf.gz"

bcftools view --threads "${THREADS}" -s "${PRO_ID},${MOM_ID},${DAD_ID}" "${MERGED_VCF}" \
| bcftools view --threads "${THREADS}" -i 'GT[0]!="0/0" && GT[0]!="./." && GT[1]=="0/0" && GT[2]=="0/0"' \
  -Oz -o "${PS_VCF}"
tabix -f -p vcf "${PS_VCF}"

bcftools view --threads "${THREADS}" -s "${PRO_ID},${MOM_ID},${DAD_ID}" "${MERGED_VCF}" \
| bcftools view --threads "${THREADS}" -i 'GT[0]!="0/0" && GT[0]!="./." && ( (GT[1]!="0/0" && GT[1]!="./.") || (GT[2]!="0/0" && GT[2]!="./.") )' \
  -Oz -o "${INH_VCF}"
tabix -f -p vcf "${INH_VCF}"

annot_one () {
  local in_vcf="$1"
  local out_prefix="$2"
  [[ -s "${in_vcf}" ]] || return 0
  "${ANNOTSV_BIN}" \
    -SVinputFile "${in_vcf}" \
    -outputFile  "${out_prefix}.tsv" \
    -genomeBuild "${GENOME_BUILD}" \
    -annotationsDir "${ANN_DIR}"
}

annot_one "${PS_VCF}"  "${TRIO_DIR}/${TRIO_TAG}.proband_specific_strict.annotsv"
annot_one "${INH_VCF}" "${TRIO_DIR}/${TRIO_TAG}.inherited_any_parent.annotsv"

echo "[DONE] Trio ${TRIO_TAG}"
echo "  Merged VCF : ${MERGED_VCF}"
echo "  Proband-only TSV : ${TRIO_DIR}/${TRIO_TAG}.proband_specific_strict.annotsv.tsv"
echo "  Inherited-any TSV: ${TRIO_DIR}/${TRIO_TAG}.inherited_any_parent.annotsv.tsv"

