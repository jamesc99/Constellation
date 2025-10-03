#!/usr/bin/env python3
import argparse
import sys
from typing import Dict, List, Optional, Tuple, Set
import pandas as pd


def normalize_gt_anyalt(gt_token: str) -> Optional[str]:
    if gt_token in (".", "./.", ".|."):
        return None
    tok = gt_token.replace("|", "/").strip()
    if tok == "":
        return None
    alleles = tok.split("/")
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    mapped: List[str] = []
    for a in alleles:
        a = a.strip()
        if a == "" or a == ".":
            return None
        if a.isdigit():
            mapped.append("0" if a == "0" else "1")
        else:
            return None
    has_ref = "0" in mapped
    has_alt = "1" in mapped
    if has_alt and not has_ref:
        return "1/1"
    if has_alt and has_ref:
        return "0/1"
    if has_ref and not has_alt:
        return "0/0"
    return None


def gt_class(gt_norm: str) -> int:
    return {"0/0": 0, "0/1": 1, "1/1": 2}[gt_norm]


def mendelian_consistent_anyalt(child_cls: int, mom_cls: int, dad_cls: int) -> bool:
    if mom_cls == 0 and dad_cls == 0:
        return child_cls == 0
    if mom_cls == 2 and dad_cls == 2:
        return child_cls == 2
    if {mom_cls, dad_cls} == {0, 2}:
        return child_cls == 1
    if mom_cls == 1 and dad_cls == 1:
        return child_cls in (0, 1, 2)
    if {mom_cls, dad_cls} == {0, 1}:
        return child_cls in (0, 1)
    if {mom_cls, dad_cls} == {1, 2}:
        return child_cls in (1, 2)
    return False


def get_gt_index(fmt: str, cache: Dict[str, int]) -> int:
    if not isinstance(fmt, str) or fmt == "":
        return -1
    if fmt in cache:
        return cache[fmt]
    parts = fmt.split(":")
    try:
        idx = parts.index("GT")
    except ValueError:
        idx = -1
    cache[fmt] = idx
    return idx


def extract_gt_from_sample(sample_field: str, gt_idx: int) -> str:
    if not isinstance(sample_field, str) or sample_field == "":
        return ""
    parts = sample_field.split(":")
    if gt_idx < 0 or gt_idx >= len(parts):
        return parts[0] if parts else ""
    return parts[gt_idx]


def autodetect_columns(df_columns: List[str]) -> Tuple[str, str, str, str]:
    cols = list(df_columns)
    idx = {c: i for i, c in enumerate(cols)}
    if all(k in idx for k in ("Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15")):
        return "Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15"
    if all(k in idx for k in ("Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12")):
        return "Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12"
    fmt_i = None
    for cand in cols:
        if "Otherinfo" in cand or cand.upper() == "FORMAT":
            fmt_i = idx[cand]
            break
    if fmt_i is None:
        fmt_i = max(0, len(cols) - 6)
    fmt_col = cols[fmt_i]
    if fmt_i + 3 < len(cols):
        p_col, m_col, f_col = cols[fmt_i + 1], cols[fmt_i + 2], cols[fmt_i + 3]
    else:
        p_col, m_col, f_col = cols[-3], cols[-2], cols[-1]
    return fmt_col, p_col, m_col, f_col


def load_gene_set(genes_file: Optional[str], genes_csv: Optional[str]) -> Optional[Set[str]]:
    gs: Set[str] = set()
    if genes_file:
        with open(genes_file) as fh:
            for line in fh:
                g = line.strip()
                if g and g != ".":
                    gs.add(g)
    if genes_csv:
        for g in genes_csv.split(","):
            g = g.strip()
            if g and g != ".":
                gs.add(g)
    return gs if gs else None


def row_overlaps_genes(gene_cell: str, gene_set: Set[str]) -> bool:
    if not isinstance(gene_cell, str):
        return False
    if gene_cell.strip() in ("", "."):
        return False
    for g in gene_cell.split(";"):
        gg = g.strip()
        if gg and gg in gene_set:
            return True
    return False


def process_trio(
    input_file: str,
    output_prefix: str,
    genes_file: Optional[str],
    genes_csv: Optional[str],
    format_col: Optional[str],
    proband_col: Optional[str],
    mother_col: Optional[str],
    father_col: Optional[str],
    af_col: str,
    af_threshold: float,
    chunksize: int,
):
    genes_set = load_gene_set(genes_file, genes_csv)
    if not genes_set:
        print("[ERROR] No genes provided. Use --genes-file or --genes.", file=sys.stderr)
        sys.exit(2)

    union_file = f"{output_prefix}_candidates_in_genes.tsv"
    denovo_file = f"{output_prefix}_denovo_in_genes.tsv"
    merr_file = f"{output_prefix}_merr_in_genes.tsv"
    all_genotyped_file = f"{output_prefix}_proband_genotyped_in_genes.tsv"

    union_header = denovo_header = merr_header = all_genotyped_header = False

    gt_index_cache: Dict[str, int] = {}
    counters = dict(total=0, af_pass=0, gene_overlap=0, proband_genotyped=0, denovo=0, merr=0)

    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)
    try:
        first_chunk = next(reader)
    except StopIteration:
        print("[INFO] Empty input.", file=sys.stderr)
        return

    if all(col is not None for col in (format_col, proband_col, mother_col, father_col)):
        fmt_c, p_c, m_c, f_c = format_col, proband_col, mother_col, father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(list(first_chunk.columns))

    def append_norm_cols(df: pd.DataFrame, trios: List[Tuple[Optional[str], Optional[str], Optional[str]]]) -> pd.DataFrame:
        s = pd.DataFrame(trios, columns=["NormGT_Proband", "NormGT_Mother", "NormGT_Father"])
        return pd.concat([df.reset_index(drop=True), s], axis=1)

    for chunk in [first_chunk] + list(reader):
        counters["total"] += len(chunk)

        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        chunk = chunk[af_vals.isna() | (af_vals < af_threshold)]
        counters["af_pass"] += int((af_vals.isna() | (af_vals < af_threshold)).sum())
        if chunk.empty:
            continue

        if "Gene.refGene" not in chunk.columns:
            continue
        gene_mask = chunk["Gene.refGene"].apply(lambda x: row_overlaps_genes(x, genes_set))
        chunk = chunk[gene_mask]
        counters["gene_overlap"] += int(gene_mask.sum())
        if chunk.empty:
            continue

        for col in (fmt_c, p_c, m_c, f_c):
            if col not in chunk.columns:
                chunk = chunk.iloc[0:0]
                break
        if chunk.empty:
            continue

        fmt_s = chunk[fmt_c].astype(str)
        p_s = chunk[p_c].astype(str)
        m_s = chunk[m_c].astype(str)
        f_s = chunk[f_c].astype(str)

        trios_norm: List[Tuple[Optional[str], Optional[str], Optional[str]]] = []
        keep_idx: List[bool] = []
        for fmt, ps, ms, fs in zip(fmt_s, p_s, m_s, f_s):
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                trios_norm.append((None, None, None))
                keep_idx.append(False)
                continue
            p_norm = normalize_gt_anyalt(extract_gt_from_sample(ps, gt_idx))
            m_norm = normalize_gt_anyalt(extract_gt_from_sample(ms, gt_idx))
            f_norm = normalize_gt_anyalt(extract_gt_from_sample(fs, gt_idx))
            keep = p_norm in ("0/0", "0/1", "1/1")
            keep_idx.append(keep)
            trios_norm.append((p_norm, m_norm, f_norm) if keep else (None, None, None))

        mask_keep = pd.Series(keep_idx, index=chunk.index)
        chunk = chunk[mask_keep]
        trios_norm = [trios_norm[i] for i, k in enumerate(keep_idx) if k]
        counters["proband_genotyped"] += int(mask_keep.sum())
        if chunk.empty:
            continue

        chunk = append_norm_cols(chunk, trios_norm)

        chunk.to_csv(all_genotyped_file, sep="\t", mode="a", header=not all_genotyped_header, index=False)
        all_genotyped_header = True

        def is_denovo(row) -> bool:
            p, m, d = row["NormGT_Proband"], row["NormGT_Mother"], row["NormGT_Father"]
            if p not in ("0/1", "1/1"):
                return False
            mom_ok = (m in ("0/0", None))
            dad_ok = (d in ("0/0", None))
            return mom_ok and dad_ok

        def is_merr(row) -> bool:
            p, m, d = row["NormGT_Proband"], row["NormGT_Mother"], row["NormGT_Father"]
            if (p not in ("0/0", "0/1", "1/1")) or (m not in ("0/0", "0/1", "1/1")) or (d not in ("0/0", "0/1", "1/1")):
                return False
            pc, mc, dc = gt_class(p), gt_class(m), gt_class(d)
            return not mendelian_consistent_anyalt(pc, mc, dc)

        denovo_df = chunk[chunk.apply(is_denovo, axis=1)].copy()
        merr_df = chunk[chunk.apply(is_merr, axis=1)].copy()

        if not denovo_df.empty:
            denovo_df.insert(0, "Flag_IsDeNovo", "1")
            denovo_df.insert(1, "Flag_IsMendelianError", "0")
            denovo_df.insert(2, "Reason", "Proband ALT; parents 0/0 or missing")
            denovo_df.to_csv(denovo_file, sep="\t", mode="a", header=not denovo_header, index=False)
            denovo_header = True
            counters["denovo"] += len(denovo_df)

        if not merr_df.empty:
            merr_df.insert(0, "Flag_IsDeNovo", "0")
            merr_df.insert(1, "Flag_IsMendelianError", "1")
            merr_df.insert(2, "Reason", "Any-ALT Mendelian inconsistency (all genotyped)")
            merr_df.to_csv(merr_file, sep="\t", mode="a", header=not merr_header, index=False)
            merr_header = True
            counters["merr"] += len(merr_df)

        if not denovo_df.empty or not merr_df.empty:
            if not denovo_df.empty and not merr_df.empty:
                merr_only = merr_df.merge(
                    denovo_df[["Chr", "Start", "End", "Ref", "Alt"]],
                    on=["Chr", "Start", "End", "Ref", "Alt"], how="left", indicator=True
                )
                merr_only = merr_only[merr_only["_merge"] == "left_only"].drop(columns=["_merge"])
                union_out = pd.concat([denovo_df, merr_only], axis=0)
            else:
                union_out = denovo_df if not denovo_df.empty else merr_df
            union_out.to_csv(union_file, sep="\t", mode="a", header=not union_header, index=False)
            union_header = True

    print("Processing complete. Summary:")
    print(f"  Total variants read:               {counters['total']}")
    print(f"  Passed AF (<{af_threshold}) filter:{counters['af_pass']}")
    print(f"  Overlapping interested genes:      {counters['gene_overlap']}")
    print(f"  Proband genotyped (kept):          {counters['proband_genotyped']} -> {all_genotyped_file}")
    print(f"  Strict de novo in genes:           {counters['denovo']} -> {denovo_file}")
    print(f"  Mendelian errors in genes:         {counters['merr']} -> {merr_file}")
    print(f"  Union candidates:                  -> {union_file}")


def main():
    ap = argparse.ArgumentParser(
        description="Find rare de novo or Mendelian-error variants in a trio ANNOVAR TSV, restricted to interested genes."
    )
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output-prefix", required=True)
    ap.add_argument("--genes-file", default=None)
    ap.add_argument("--genes", default=None)
    ap.add_argument("--format-col", default=None)
    ap.add_argument("--proband-col", default=None)
    ap.add_argument("--mother-col", default=None)
    ap.add_argument("--father-col", default=None)
    ap.add_argument("--af-col", default="gnomad41_genome_AF")
    ap.add_argument("--af-threshold", type=float, default=0.001)
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()

    process_trio(
        input_file=args.input,
        output_prefix=args.output_prefix,
        genes_file=args.genes_file,
        genes_csv=args.genes,
        format_col=args.format_col,
        proband_col=args.proband_col,
        mother_col=args.mother_col,
        father_col=args.father_col,
        af_col=args.af_col,
        af_threshold=args.af_threshold,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()

