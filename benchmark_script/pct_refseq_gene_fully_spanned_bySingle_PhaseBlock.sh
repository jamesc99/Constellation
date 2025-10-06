#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <whatshap_block_list> <refseq_gene_bed> [output_file]"
  exit 1
fi

BLOCK_LIST="$1"
GENES_BED="$2"
OUT="${3:-pct_fully_spanned_genes_coveredBySingle_phasedBlock.txt}"

if ! command -v bedtools >/dev/null 2>&1; then
  echo "bedtools not found in PATH"; exit 1
fi

FULLY_PHASED_GENES=$(
  awk 'NR==1{
         for(i=1;i<=NF;i++){
           if($i=="chrom"||$i=="#chrom"||$i=="CHR"){c=i}
           if($i=="from"||$i=="start"||$i=="FROM"){f=i}
           if($i=="to"||$i=="end"||$i=="TO"){t=i}
         }
         next
       }
       {print $c "\t" $f "\t" $t}' "$BLOCK_LIST" \
  | bedtools intersect -a "$GENES_BED" -b stdin -f 1.0 -u \
  | wc -l
)

TOTAL_GENES=$(wc -l < "$GENES_BED")
if [[ "${TOTAL_GENES}" -eq 0 ]]; then
  PCT="0"
else
  PCT=$(echo "scale=4; $FULLY_PHASED_GENES / $TOTAL_GENES * 100" | bc)
fi

{
  echo "Total number of RefSeq genes: $TOTAL_GENES"
  echo "Number of genes fully spanned by a single phase block: $FULLY_PHASED_GENES"
  echo "Percentage of spanned genes: $PCT %"
} > "$OUT"

echo "Successfully wrote statistics to $OUT"

