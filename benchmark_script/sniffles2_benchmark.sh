#!/bin/bash
#SBATCH --job-name=sniffles2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

INPUT_BAM="${1:?Usage: $0 <input.bam>}"
THREADS="${THREADS:-8}"
REF_FASTA="${REF_FASTA:?Set REF_FASTA to reference FASTA path}"
SNIFFLES_BIN="${SNIFFLES_BIN:-sniffles}"

time_to_seconds() {
  local t="$1"; local m=0; local s=0
  if [[ "$t" == *m* ]]; then m="${t%m*}"; s="${t#*m}"; s="${s%s}"; else s="${t%s}"; fi
  echo "$m * 60 + $s" | bc
}

NAME="${SLURM_JOB_ID:-$(basename "${INPUT_BAM%%.bam}")}"
TIME_LOG_TMP="$(mktemp)"

{ time "${SNIFFLES_BIN}" --input "${INPUT_BAM}" --vcf "${NAME}.vcf" --sample-id "${NAME}" --reference "${REF_FASTA}" --threads "${THREADS}" > "${INPUT_BAM%.bam}.log" 2>&1; } 2> "${TIME_LOG_TMP}"

REAL_TIME="$(grep '^real' "${TIME_LOG_TMP}" | awk '{print $2}')"
USER_TIME_STR="$(grep '^user' "${TIME_LOG_TMP}" | awk '{print $2}')"
SYS_TIME_STR="$(grep '^sys' "${TIME_LOG_TMP}" | awk '{print $2}')"
USER_SEC="$(time_to_seconds "${USER_TIME_STR}")"
SYS_SEC="$(time_to_seconds "${SYS_TIME_STR}")"
TOTAL_CPU_SEC="$(echo "${USER_SEC} + ${SYS_SEC}" | bc)"

{
  echo "job_id: ${NAME}"
  echo "total_cpu_time: ${TOTAL_CPU_SEC}s (${USER_SEC}s user + ${SYS_SEC}s sys)"
  echo "wall_clock_time: ${REAL_TIME}"
  echo "---"
} >> running_time.log

rm -f "${TIME_LOG_TMP}"

