#!/bin/bash

# fastp_array.sh — Run fastp on paired-end FASTQs listed as R1/R2 on alternating lines.
#
# Usage:
#   N=$(( $(wc -l < fastqlist.tsv) / 2 ))
#   sbatch --array=1-"$N" fastp_array.sh fastqlist.tsv trimmed
#
#SBATCH --job-name=fastp
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=mhgcp
#SBATCH --time=40:00:00
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

set -euo pipefail
IFS=$'\n\t'

LIST_FILE=${1:?Provide fastqlist.tsv}
OUTDIR=${2:-trimmed}

# ---------------- 1) Robust TMPDIR on scratch (prevents "Unable to create TMPDIR") ----------------
# Pick a real, big scratch location on your cluster. Adjust if needed:
SCRATCH_ROOT="/research/sedlazeck/ryan/tmp"
export TMPDIR="${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"

# ---------------- 2) Locate fastp ----------------
# Preferred: absolute path (most robust on non-interactive Slurm jobs)
FASTP_BIN="/research/sedlazeck/ryan/miniconda3/bin/fastp"

# Fallback: if you prefer conda activation, uncomment this block and set ENV_NAME properly.
# CONDA_ROOT="/research/sedlazeck/ryan/miniconda3"
# ENV_NAME="base"   # or a specific env that has fastp installed
# if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
#   source "${CONDA_ROOT}/etc/profile.d/conda.sh"
#   conda activate "${ENV_NAME}"
#   FASTP_BIN="fastp"
# else
#   export PATH="${CONDA_ROOT}/bin:${PATH}"
# fi

command -v "$FASTP_BIN" >/dev/null 2>&1 || { echo "ERROR: fastp not found at $FASTP_BIN"; exit 1; }

# ---------------- 3) Read the current array pair (R1/R2 on alternating lines) ----------------
total_lines=$(wc -l < "$LIST_FILE")
if (( total_lines % 2 != 0 )); then
  echo "ERROR: $LIST_FILE has an odd number of lines ($total_lines). Each pair must be R1 then R2."
  exit 1
fi

PAIR_INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}
start_line=$(( (PAIR_INDEX - 1) * 2 + 1 ))
r1=$(sed -n "${start_line}p" "$LIST_FILE")
r2=$(sed -n "$(( start_line + 1 ))p" "$LIST_FILE")

if [[ -z "${r1}" || -z "${r2}" ]]; then
  echo "ERROR: Could not read lines ${start_line} and $((start_line+1)) from $LIST_FILE"
  exit 1
fi
[[ -f "$r1" ]] || { echo "ERROR: R1 missing: $r1"; exit 1; }
[[ -f "$r2" ]] || { echo "ERROR: R2 missing: $r2"; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/reports"

# ---------------- 4) Derive sample name from R1 ----------------
base1=$(basename "$r1")
name="${base1}"
# strip extensions
name="${name%.fastq.gz}"; name="${name%.fq.gz}"; name="${name%.fastq}"; name="${name%.fq}"
# strip common R1 tokens
if [[ "$name" == *"_R1_"* ]]; then
  name="${name%%_R1_*}"
elif [[ "$name" == *"_R1"* ]]; then
  name="${name%%_R1*}"
elif [[ "$name" == *".R1."* ]]; then
  name="${name%%.R1.*}"
elif [[ "$name" == *"_1_"* ]]; then
  name="${name%%_1_*}"
elif [[ "$name" == *"_1"* ]]; then
  name="${name%%_1*}"
fi

out_r1="$OUTDIR/${name}_R1.trimmed.fastq.gz"
out_r2="$OUTDIR/${name}_R2.trimmed.fastq.gz"
json="$OUTDIR/reports/${name}.json"
html="$OUTDIR/reports/${name}.html"

threads=${SLURM_CPUS_PER_TASK:-8}

echo "[$(date)] Pair #$PAIR_INDEX"
echo "R1: $r1"
echo "R2: $r2"
echo "Sample: $name"
echo "Threads: $threads"
echo "TMPDIR: $TMPDIR"
echo

# ---------------- 5) Run fastp ----------------
# Adapter sequences: Illumina universal (R1) and reverse-complement for R2 commonly used.
"$FASTP_BIN" \
  --dedup --dup_calc_accuracy 5 \
  -i "$r1" -I "$r2" \
  -o "$out_r1" -O "$out_r2" \
  --thread "$threads" \
  --json "$json" \
  --html "$html"

echo "[$(date)] Done: $name"

