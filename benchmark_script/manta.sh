#!/bin/bash
#SBATCH --job-name=manta_single
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

BAM_OR_CRAM="${1:?Usage: $0 <input.bam|input.cram>}"
MANTA_INSTALL_DIR="${MANTA_INSTALL_DIR:?Set MANTA_INSTALL_DIR (path to manta install, e.g. /opt/manta-1.6.0)}"
REF_FASTA="${REF_FASTA:?Set REF_FASTA to the reference FASTA path (matching the BAM/CRAM)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

SAMPLE_NAME="$(basename "${BAM_OR_CRAM}")"
SAMPLE_NAME="${SAMPLE_NAME%.bam}"
SAMPLE_NAME="${SAMPLE_NAME%.cram}"

OUTDIR="${SAMPLE_NAME}_manta"
mkdir -p "${OUTDIR}"

"${MANTA_INSTALL_DIR}/bin/configManta.py" \
  --bam "${BAM_OR_CRAM}" \
  --referenceFasta "${REF_FASTA}" \
  --runDir "${OUTDIR}"

cd "${OUTDIR}"
./runWorkflow.py -j "${THREADS}"

VCF="results/variants/candidateSV.vcf.gz"
[[ -f "${VCF}" ]] || { echo "ERROR: Expected VCF not found: ${VCF}"; exit 1; }
echo "Done. VCF: $(realpath "${VCF}")"

