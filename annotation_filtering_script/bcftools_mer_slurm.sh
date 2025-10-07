BATCH -J mendel_MER
#SBATCH -o mendel_MER_%A_%a.out
#SBATCH -e mendel_MER_%A_%a.err
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 3:00:00
#SBATCH -p compute

set -euo pipefail

TSV="${1:?Usage: sbatch --array=1-N bcftools_mer_slurm.sh <trios.tsv> [outdir]}"
OUTDIR="${2:-mer_results}"
MER_SCRIPT="${MER_SCRIPT:?Set MER_SCRIPT to a script that converts mendelian2 summary -> MER TSV}"

mkdir -p "${OUTDIR}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
LINE_NUM="${SLURM_ARRAY_TASK_ID:?Set --array=1-N when submitting}"
LINE="$(sed -n "${LINE_NUM}p" "${TSV}")"
[[ -n "${LINE}" ]] || { echo "No line ${LINE_NUM} in ${TSV}" >&2; exit 2; }

PROBAND_VCF="$(awk '{print $1}' <<<"${LINE}")"
MOTHER_VCF="$(awk '{print $2}' <<<"${LINE}")"
FATHER_VCF="$(awk '{print $3}' <<<"${LINE}")"
SEX="$(awk '{print $4}' <<<"${LINE}")"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
TABIX_BIN="${TABIX_BIN:-tabix}"

for f in "${PROBAND_VCF}" "${MOTHER_VCF}" "${FATHER_VCF}"; do
  [[ -s "${f}" ]] || { echo "Missing VCF: ${f}" >&2; exit 3; }
  [[ -s "${f}.tbi" || -s "${f}.csi" ]] || { echo "Missing index for ${f}" >&2; exit 4; }
done
[[ "${SEX}" == "1" || "${SEX}" == "2" ]] || { echo "SEX must be 1 or 2, got: ${SEX}" >&2; exit 5; }

PROBAND_ID="$("${BCFTOOLS_BIN}" query -l "${PROBAND_VCF}" | head -n1 || true)"
MOTHER_ID="$("${BCFTOOLS_BIN}" query -l "${MOTHER_VCF}" | head -n1 || true)"
FATHER_ID="$("${BCFTOOLS_BIN}" query -l "${FATHER_VCF}" | head -n1 || true)"

if [[ -n "${PROBAND_ID}" ]]; then
  PREFIX="${OUTDIR}/${PROBAND_ID}"
else
  BASE="$(basename "${PROBAND_VCF}")"
  PREFIX="${OUTDIR}/${BASE%.vcf.gz}"
fi

SEX_PREFIX=$([[ "${SEX}" == "1" ]] && echo "1X:" || echo "2X:")

"${BCFTOOLS_BIN}" merge --threads "${THREADS}" -m all -Oz \
  -o "${PREFIX}.merged.vcf.gz" \
  "${PROBAND_VCF}" "${MOTHER_VCF}" "${FATHER_VCF}"

"${TABIX_BIN}" -p vcf "${PREFIX}.merged.vcf.gz"

"${BCFTOOLS_BIN}" +mendelian2 "${PREFIX}.merged.vcf.gz" \
  --pfm "${SEX_PREFIX}${PROBA_ID},${FATHER_ID},${MOTHER_ID}" \
  --mode c > "${PREFIX}.summary.txt"

# --- compute MER from the summary file (your script) ---
bash "$MER_SCRIPT" "${PREFIX}.summary.txt" > "${PREFIX}.mer.tsv"

echo "Done."
echo "  Summary: ${PREFIX}.summary.txt"
echo "  MER:     ${PREFIX}.mer.tsv"

