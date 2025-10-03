BATCH --job-name=whatshap_stats
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=3:00:00

set -euo pipefail

INPUT_VCFGZ="${1:?Usage: $0 <input.vcf.gz> <name>}"
NAME="${2:?Usage: $0 <input.vcf.gz> <name>}"

CHR_LENGTHS="${CHR_LENGTHS:?Set CHR_LENGTHS to a chr_lengths.txt file}"
BLOCK_LIST="${BLOCK_LIST:-whatshap_stat_phased_blocks.bed}"

whatshap stats \
  --block-list="${BLOCK_LIST}" \
  --tsv="${NAME}_stats.tsv" \
  --chr-lengths="${CHR_LENGTHS}" \
  "${INPUT_VCFGZ}"

