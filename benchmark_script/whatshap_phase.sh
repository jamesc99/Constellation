BATCH --job-name=whatshap_phasing
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

INPUT_VCF="${1:?Usage: $0 <input.vcf.gz> <input.bam>}"
INPUT_BAM="${2:?Usage: $0 <input.vcf.gz> <input.bam>}"
REFERENCE="${REFERENCE:?Set REFERENCE to reference FASTA path}"
OUTPUT_DIR="${OUTPUT_DIR:-$(pwd)}"

whatshap phase \
  --reference "${REFERENCE}" \
  -o "${OUTPUT_DIR}/phased.vcf.gz" \
  --ignore-read-groups \
  "${INPUT_VCF}" "${INPUT_BAM}"

