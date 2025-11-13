#SBATCH --job-name=clair3
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

INPUT_BAM="${1:?Usage: $0 <input.bam>}"

# --- Config (override with environment variables) ---
REFERENCE="${REFERENCE:?Set REFERENCE to GRCh38 FASTA path}"
CLAIR3_MODEL="${CLAIR3_MODEL:?Set CLAIR3_MODEL to model directory path}"
THREADS="${THREADS:-12}"
OUTPUT_DIR="${OUTPUT_DIR:-$(pwd)}"
PLATFORM="${PLATFORM:-ont}"

run_clair3.sh \
  --bam_fn="${INPUT_BAM}" \
  --ref_fn="${REFERENCE}" \
  --threads="${THREADS}" \
  --platform="${PLATFORM}" \
  --model_path="${CLAIR3_MODEL}" \
  --output="${OUTPUT_DIR}"
