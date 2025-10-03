#!/bin/bash
#SBATCH --job-name=rtg_vcfeval
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=3:00:00
#SBATCH --partition=compute

set -euo pipefail

INPUT_VCFGZ="${1:?Usage: $0 <comp.vcf.gz> <name>}"
NAME="${2:?Usage: $0 <comp.vcf.gz> <name>}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
RTG_BIN="${RTG_BIN:-rtg}"
THREADS="${THREADS:-4}"
REF_SDF="${REF_SDF:?Set REF_SDF to the path of the reference .sdf directory}"

T2T_SNV_VCF="${T2T_SNV_VCF:?Set T2T_SNV_VCF to T2T-Q100 v1.1 SNV-only VCF path}"
T2T_SNVINDEL_VCF="${T2T_SNVINDEL_VCF:?Set T2T_SNVINDEL_VCF to T2T-Q100 v1.1 SNVs+indels VCF path}"
T2T_BED="${T2T_BED:?Set T2T_BED to T2T-Q100 v1.1 benchmark BED path}"

NIST421_SNV_VCF="${NIST421_SNV_VCF:?Set NIST421_SNV_VCF to NIST v4.2.1 SNV-only VCF path}"
NIST421_SNVINDEL_VCF="${NIST421_SNVINDEL_VCF:?Set NIST421_SNVINDEL_VCF to NIST v4.2.1 SNVs+indels VCF path}"
NIST421_BED="${NIST421_BED:?Set NIST421_BED to NIST v4.2.1 benchmark BED path}"

CMRG_SNV_VCF="${CMRG_SNV_VCF:?Set CMRG_SNV_VCF to CMRG v1.0 SNV-only VCF path}"
CMRG_SNVINDEL_VCF="${CMRG_SNVINDEL_VCF:?Set CMRG_SNVINDEL_VCF to CMRG v1.0 SNVs+indels VCF path}"
CMRG_BED="${CMRG_BED:?Set CMRG_BED to CMRG v1.0 benchmark BED path}"

"${BCFTOOLS_BIN}" view --threads "${THREADS}" -v snps "${INPUT_VCFGZ}" -Oz -o "${NAME}_snvs_only.vcf.gz"
"${BCFTOOLS_BIN}" index -ft --threads "${THREADS}" "${NAME}_snvs_only.vcf.gz"

mkdir -p "${NAME}_rtg_output"

run_vcfeval() {
  BASE_VCF="$1"
  BED_FILE="$2"
  COMP_VCF="$3"
  OUT_PREFIX="$4"

  "${RTG_BIN}" vcfeval --threads="${THREADS}" --ref-overlap \
    -b "${BASE_VCF}" \
    --bed-regions "${BED_FILE}" \
    -c "${COMP_VCF}" \
    -t "${REF_SDF}" \
    -o "./${NAME}_rtg_output/${OUT_PREFIX}_GTmatch" \
    --vcf-score-field QUAL

  "${RTG_BIN}" vcfeval --threads="${THREADS}" --ref-overlap --squash-ploidy \
    -b "${BASE_VCF}" \
    --bed-regions "${BED_FILE}" \
    -c "${COMP_VCF}" \
    -t "${REF_SDF}" \
    -o "./${NAME}_rtg_output/${OUT_PREFIX}_NoGT" \
    --vcf-score-field QUAL
}

run_vcfeval "${T2T_SNV_VCF}" "${T2T_BED}" "${NAME}_snvs_only.vcf.gz" "nistq100_snv_only_rtg_vcfeval_output"
run_vcfeval "${T2T_SNVINDEL_VCF}" "${T2T_BED}" "${INPUT_VCFGZ}" "nistq100_snv_indels_rtg_vcfeval_output"

run_vcfeval "${NIST421_SNV_VCF}" "${NIST421_BED}" "${NAME}_snvs_only.vcf.gz" "nistv4.2.1_snv_only_rtg_vcfeval_output"
run_vcfeval "${NIST421_SNVINDEL_VCF}" "${NIST421_BED}" "${INPUT_VCFGZ}" "nistv4.2.1_snv_indels_rtg_vcfeval_output"

run_vcfeval "${CMRG_SNV_VCF}" "${CMRG_BED}" "${NAME}_snvs_only.vcf.gz" "cmrgv1.0.0_snv_only_rtg_vcfeval_output"
run_vcfeval "${CMRG_SNVINDEL_VCF}" "${CMRG_BED}" "${INPUT_VCFGZ}" "cmrgv1.0.0_snv_indels_rtg_vcfeval_output"

