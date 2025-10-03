#!/usr/bin/env bash
set -euo pipefail

stats_file="${1:-stats.txt}"

# Header
echo -e "trio_id\tchild\tfather\tmother\tngood\tnmerr\tnmissing\tMER(%)\tMissingRate(%)"

awk '
BEGIN{FS="\t"; OFS="\t"}
/^TRIO[[:space:]]+/{
  id=$2; child=$3; father=$4; mother=$5
  ch[id]=child; fa[id]=father; mo[id]=mother; ids[id]=1
}
/^ngood[[:space:]]+/{
  for(i=2;i<=NF;i++) ng[i-1]=$i
}
/^nmerr[[:space:]]+/{
  for(i=2;i<=NF;i++) me[i-1]=$i
}
/^nmissing[[:space:]]+/{
  for(i=2;i<=NF;i++) mi[i-1]=$i
}
END{
  t_ng=0; t_me=0; t_mi=0
  # ids[] may not be set if TRIO lines are absent; infer count from ng[]
  # Build a numeric index list 1..N based on ng[]
  n=0
  for (k in ng){ if (k>n) n=k }
  for (i=1;i<=n;i++){
    ngood = (i in ng)? ng[i] : 0
    nmerr = (i in me)? me[i] : 0
    nmiss = (i in mi)? mi[i] : 0
    denom = ngood + nmerr
    total = denom + nmiss
    mer   = (denom>0)? (100.0*nmerr/denom) : 0.0
    missr = (total>0)? (100.0*nmiss/total) : 0.0
    t_ng += ngood; t_me += nmerr; t_mi += nmiss
    cid = (i in ch)? ch[i] : "-"
    cfa = (i in fa)? fa[i] : "-"
    cmo = (i in mo)? mo[i] : "-"
    printf "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\n", i, cid, cfa, cmo, ngood, nmerr, nmiss, mer, missr
  }
  t_denom = t_ng + t_me
  t_total = t_denom + t_mi
  t_mer   = (t_denom>0)? (100.0*t_me/t_denom) : 0.0
  t_missr = (t_total>0)? (100.0*t_mi/t_total) : 0.0
  printf "TOTAL\t-\t-\t-\t%d\t%d\t%d\t%.2f\t%.2f\n", t_ng, t_me, t_mi, t_mer, t_missr
}' "$stats_file"

