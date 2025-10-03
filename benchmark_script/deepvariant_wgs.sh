BATCH --job-name=deepvariant_wgs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=12
#SBATCH --mem=80G

set -euo pipefail

BAM_IN="${1:?Usage: $0 <sample.bam>}"
REF_FASTA="${REF_FASTA:?Set REF_FASTA to reference FASTA path}"
BIN_VERSION="${BIN_VERSION:-1.9.0}"
OUT_DIR="${OUT_DIR:-${PWD}/output}"
INT_DIR="${INT_DIR:-${OUT_DIR}/intermediate_results_dir}"
APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${PWD}/.apptainer_cache}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-12}}"

mkdir -p "${OUT_DIR}" "${INT_DIR}" "${APPTAINER_CACHEDIR}"

BAM_DIR="$(dirname "${BAM_IN}")"
BAM_NAME="$(basename "${BAM_IN}")"
SAMPLE="${BAM_NAME%.bam}"
REF_DIR="$(dirname "${REF_FASTA}")"
REF_BASE="$(basename "${REF_FASTA}")"
JOB_TMP="${JOB_TMP:-/tmp/deepvariant.${SLURM_JOB_ID:-$$}}"
mkdir -p "${JOB_TMP}"

[[ -f "${BAM_IN}" ]] || { echo "[ERROR] BAM not found: ${BAM_IN}"; exit 1; }
[[ -f "${BAM_IN}.bai" ]] || { echo "[ERROR] BAM index not found: ${BAM_IN}.bai"; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "[ERROR] Reference FASTA not found: ${REF_FASTA}"; exit 1; }
[[ -f "${REF_FASTA}.fai" ]] || { echo "[ERROR] Missing FASTA index: ${REF_FASTA}.fai"; exit 1; }

export APPTAINER_CACHEDIR
apptainer run \
  -B "${JOB_TMP}:${JOB_TMP}","${BAM_DIR}":/mnt/input,"${OUT_DIR}":/mnt/output,"${REF_DIR}":/mnt/ref \
  docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type WGS \
    --ref /mnt/ref/"${REF_BASE}" \
    --reads /mnt/input/"${BAM_NAME}" \
    --output_vcf /mnt/output/"${SAMPLE}.deepvariant.vcf.gz" \
    --output_gvcf /mnt/output/"${SAMPLE}.deepvariant.g.vcf.gz" \
    --intermediate_results_dir /mnt/output/intermediate_results_dir \
    --num_shards "${THREADS}"

echo "[DONE] ${OUT_DIR}/${SAMPLE}.deepvariant.vcf.gz"
echo "[DONE] ${OUT_DIR}/${SAMPLE}.deepvariant.g.vcf.gz"

