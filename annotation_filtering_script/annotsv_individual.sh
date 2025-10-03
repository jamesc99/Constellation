#!/usr/bin/env bash
#SBATCH --job-name=annotsv_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --array=1-2

set -euo pipefail

ANNOTSV_DIR="${ANNOTSV_DIR:?Set ANNOTSV_DIR to the AnnotSV install dir}"
ANNOTSV_BIN="${ANNOTSV_BIN:-${ANNOTSV_DIR}/bin/AnnotSV}"
ANN_DIR="${ANN_DIR:?Set ANN_DIR to the AnnotSV annotations directory}"

LIST_FILE="${LIST_FILE:-vcflist.txt}"
OUT_DIR="${OUT_DIR:-annotsv_out}"
EMIT_VCF="${EMIT_VCF:-0}"
GENES_FILE="${GENES_FILE:-}"

[[ -f "${LIST_FILE}" ]] || { echo "ERROR: LIST_FILE not found: ${LIST_FILE}" >&2; exit 1; }
N_LINES="$(wc -l < "${LIST_FILE}")"
[[ "${SLURM_ARRAY_TASK_ID:-0}" -ge 1 && "${SLURM_ARRAY_TASK_ID}" -le "${N_LINES}" ]] || {
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-unset} out of range (1..${N_LINES})" >&2; exit 1; }

IN_VCF="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST_FILE}" | tr -d '\r')"
[[ -n "${IN_VCF}" && -f "${IN_VCF}" ]] || { echo "ERROR: Input VCF missing on line ${SLURM_ARRAY_TASK_ID}: '${IN_VCF}'" >&2; exit 1; }

base="$(basename "${IN_VCF}")"
base="${base%.vcf.gz}"; base="${base%.vcf}"; base="${base%.gz}"; base="${base%.sv}"
if [[ "$(awk -F'_' '{print NF}' <<<"${base}")" -ge 2 ]]; then
  OUT_PREFIX="$(awk -F'_' '{print $1"_"$2}' <<<"${base}")"
else
  OUT_PREFIX="${base}"
fi

mkdir -p "${OUT_DIR}"
SUBDIR="${OUT_DIR}/${OUT_PREFIX}"
mkdir -p "${SUBDIR}"
export TMPDIR="${SUBDIR}/tmp"
mkdir -p "${TMPDIR}"

OUT_BASENAME="${SUBDIR}/${OUT_PREFIX}.annotsv"

[[ -x "${ANNOTSV_BIN}" ]] || { echo "ERROR: AnnotSV not found at ${ANNOTSV_BIN}" >&2; exit 1; }
[[ -d "${ANN_DIR}" ]] || { echo "ERROR: annotations dir not found: ${ANN_DIR}" >&2; exit 1; }
[[ -z "${GENES_FILE}" || -f "${GENES_FILE}" ]] || { echo "ERROR: GENES_FILE not found: ${GENES_FILE}" >&2; exit 1; }

CMD=( "${ANNOTSV_BIN}"
  -SVinputFile "${IN_VCF}"
  -outputFile "${OUT_BASENAME}.tsv"
  -genomeBuild GRCh38
  -annotationsDir "${ANN_DIR}"
)
[[ -n "${GENES_FILE}" ]] && CMD+=( -candidateGenesFile "${GENES_FILE}" -candidateGenesFileAnnotations 1 )

echo "[INFO] [${OUT_PREFIX}] Running AnnotSV"
"${CMD[@]}"

if [[ "${EMIT_VCF}" == "1" ]]; then
  "${ANNOTSV_BIN}" \
    -SVinputFile "${IN_VCF}" \
    -outputFile "${OUT_BASENAME}.vcf" \
    -genomeBuild GRCh38 \
    -annotationsDir "${ANN_DIR}" \
    -outputAnnotations 1
  if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    bgzip -f "${OUT_BASENAME}.vcf"
    tabix -f -p vcf "${OUT_BASENAME}.vcf.gz"
  fi
fi

touch "${SUBDIR}/_DONE"

echo "[DONE] [${OUT_PREFIX}]"
echo "  TSV : ${OUT_BASENAME}.tsv"
[[ "${EMIT_VCF}" == "1" && -f "${OUT_BASENAME}.vcf.gz" ]] && echo "  VCF : ${OUT_BASENAME}.vcf.gz"
echo "  TMP : ${TMPDIR}"

