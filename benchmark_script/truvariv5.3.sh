BATCH --job-name=truvari53
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

INPUT_VCFGZ="${1:?Usage: $0 <comp.vcf.gz> <name>}"
NAME="${2:?Usage: $0 <comp.vcf.gz> <name>}"

# --- Config (override with environment variables) ---
TRUVARI_BIN="${TRUVARI_BIN:-truvari}"
REF_FASTA="${REF_FASTA:?Set REF_FASTA to GRCh38 FASTA path}"
T2T_BASE_VCF="${T2T_BASE_VCF:?Set T2T_BASE_VCF to the GIAB T2T-Q100 v1.1 SV VCF path}"
T2T_INCLUDE_BED="${T2T_INCLUDE_BED:?Set T2T_INCLUDE_BED to the GIAB T2T-Q100 v1.1 benchmark BED path}"
CMRG_BASE_VCF="${CMRG_BASE_VCF:?Set CMRG_BASE_VCF to the CMRG v1.0 SV VCF path}"
CMRG_INCLUDE_BED="${CMRG_INCLUDE_BED:?Set CMRG_INCLUDE_BED to the CMRG v1.0 benchmark BED path}"

THREADS="${THREADS:-8}"
REFDIST="${REFDIST:-2000}"
PCTSEQ="${PCTSEQ:-0.7}"
PCTSIZE="${PCTSIZE:-0.7}"
PCTOVL="${PCTOVL:-0.0}"
SIZEMIN="${SIZEMIN:-50}"
SIZEFILT="${SIZEFILT:-50}"
SIZEMAX="${SIZEMAX:-50000}"
PICK="${PICK:-ac}"
EXTEND="${EXTEND:-0}"
CHUNKSIZE="${CHUNKSIZE:-5000}"
ALIGNER="${ALIGNER:-mafft}"

"${TRUVARI_BIN}" version

echo "GIAB T2T_Q100_v1.1 results"
"${TRUVARI_BIN}" bench \
  --base "${T2T_BASE_VCF}" \
  --comp "${INPUT_VCFGZ}" \
  --output "./${NAME}_q100_truvari_v5.3/" \
  --includebed "${T2T_INCLUDE_BED}" \
  --refdist "${REFDIST}" \
  --pctseq "${PCTSEQ}" \
  --pctsize "${PCTSIZE}" \
  --pctovl "${PCTOVL}" \
  --passonly \
  --sizemin "${SIZEMIN}" \
  --sizefilt "${SIZEFILT}" \
  --sizemax "${SIZEMAX}" \
  --pick "${PICK}" \
  --extend "${EXTEND}" \
  --chunksize "${CHUNKSIZE}"

"${TRUVARI_BIN}" refine \
  --reference "${REF_FASTA}" \
  --regions "./${NAME}_q100_truvari_v5.3/candidate.refine.bed" \
  --use-original-vcfs \
  --coords R \
  --align "${ALIGNER}" \
  --threads "${THREADS}" \
  "./${NAME}_q100_truvari_v5.3/"

echo "CMRG v1.0 results"
"${TRUVARI_BIN}" bench \
  --base "${CMRG_BASE_VCF}" \
  --comp "${INPUT_VCFGZ}" \
  --output "./${NAME}_cmrg_truvari_v5.3/" \
  --includebed "${CMRG_INCLUDE_BED}" \
  --refdist "${REFDIST}" \
  --pctseq "${PCTSEQ}" \
  --pctsize "${PCTSIZE}" \
  --pctovl "${PCTOVL}" \
  --passonly \
  --sizemin "${SIZEMIN}" \
  --sizefilt "${SIZEFILT}" \
  --sizemax "${SIZEMAX}" \
  --pick "${PICK}" \
  --extend "${EXTEND}" \
  --chunksize "${CHUNKSIZE}"

"${TRUVARI_BIN}" refine \
  --reference "${REF_FASTA}" \
  --regions "./${NAME}_cmrg_truvari_v5.3/candidate.refine.bed" \
  --use-original-vcfs \
  --coords R \
  --align "${ALIGNER}" \
  --threads "${THREADS}" \
  "./${NAME}_cmrg_truvari_v5.3/"

