#!/bin/bash
#SBATCH --job-name=bwa_mapping
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --partition=medium

set -euo pipefail

READ1="${1:?Usage: $0 <read1.fq.gz> <read2.fq.gz> <sample_name> [platform]}"
READ2="${2:?Usage: $0 <read1.fq.gz> <read2.fq.gz> <sample_name> [platform]}"
SAMPLE_NAME="${3:?Usage: $0 <read1.fq.gz> <read2.fq.gz> <sample_name> [platform]}"
PLATFORM_NAME="${4:-Illumina}"

BWA_BIN="${BWA_BIN:-bwa}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
REF_FASTA="${REF_FASTA:?Set REF_FASTA to reference FASTA path}"

THREADS="${THREADS:-16}"

"${BWA_BIN}" mem -t "${THREADS}" -M \
  -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM_NAME}" \
  "${REF_FASTA}" "${READ1}" "${READ2}" > "${SAMPLE_NAME}.sam"

"${SAMTOOLS_BIN}" view -u -@"${THREADS}" -o "${SAMPLE_NAME}.bam" "${SAMPLE_NAME}.sam"
"${SAMTOOLS_BIN}" sort -@"${THREADS}" -o "${SAMPLE_NAME}_sorted.bam" "${SAMPLE_NAME}.bam"
"${SAMTOOLS_BIN}" index -@"${THREADS}" "${SAMPLE_NAME}_sorted.bam"

rm -f "${SAMPLE_NAME}.sam" "${SAMPLE_NAME}.bam"

