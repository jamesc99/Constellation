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


def parse_gene_cell(gene_cell: str) -> List[str]:
    if not isinstance(gene_cell, str) or gene_cell.strip() in ("", "."):
        return []
    return [g.strip() for g in gene_cell.split(";") if g.strip() not in ("", ".")]


def load_gene_set(gene_list_arg: Optional[str]) -> Set[str]:
    if gene_list_arg is None:
        return {"DMD"}
    try:
        with open(gene_list_arg, "r") as fh:
            genes = {line.strip().upper() for line in fh if line.strip() and not line.startswith("#")}
            if genes:
                return genes
    except Exception:
        pass
    parts = [p.strip().upper() for p in gene_list_arg.split(",") if p.strip()]
    return set(parts)


def is_mendelian_error(p_norm: Optional[str], m_norm: Optional[str], f_norm: Optional[str]) -> bool:
    if p_norm is None or m_norm is None or f_norm is None:
        return False
    if (m_norm == "0/0" and f_norm == "0/0") and (p_norm in ("0/1", "1/1")):
        return True
    if (m_norm == "1/1" and f_norm == "1/1") and (p_norm == "0/0"):
        return True
    if p_norm == "1/1" and not (m_norm in ("0/1", "1/1") and f_norm in ("0/1", "1/1")):
        return True
    return False


def autodetect_columns(df_columns: List[str]) -> Tuple[str, str, str, str]:
    cols = list(df_columns)
    name_to_idx = {c: i for i, c in enumerate(cols)}
    if all(k in name_to_idx for k in ("Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15")):
        return "Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15"
    if all(k in name_to_idx for k in ("Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12")):
        return "Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12"
    fmt_idx = None
    for cand in cols:
        if "Otherinfo" in cand or cand.upper() == "FORMAT":
            fmt_idx = name_to_idx[cand]
            break
    if fmt_idx is None:
        fmt_idx = max(0, len(cols) - 6)
    fmt_col = cols[fmt_idx]
    if fmt_idx + 3 < len(cols):
        p_col, m_col, f_col = cols[fmt_idx + 1], cols[fmt_idx + 2], cols[fmt_idx + 3]
    else:
        p_col, m_col, f_col = cols[-3], cols[-2], cols[-1]
    return fmt_col, p_col, m_col, f_col


def process_trio(
    input_file: str,
    output_prefix: str,
    genes_of_interest: Set[str],
    format_col: Optional[str] = None,
    proband_col: Optional[str] = None,
    mother_col: Optional[str] = None,
    father_col: Optional[str] = None,
    af_col: str = "gnomad41_genome_AF",
    af_threshold: float = 0.001,
    allow_missing_af: bool = True,
    chunksize: int = 100000,
):
    denovo_out = f"{output_prefix}_denovo.tsv"
    comphet_out = f"{output_prefix}_comphet.tsv"

    header_cols: List[str] = []
    gt_index_cache: Dict[str, int] = {}
    denovo_rows: List[Dict[str, str]] = []
    gene_maternal: Dict[str, Set[Tuple]] = {}
    gene_paternal: Dict[str, Set[Tuple]] = {}
    row_store: Dict[Tuple, Dict[str, str]] = {}
    row_order: Dict[Tuple, int] = {}
    order_counter = 0

    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)
    try:
        first_chunk = next(reader)
    except StopIteration:
        pd.DataFrame(columns=[]).to_csv(denovo_out, sep="\t", index=False)
        pd.DataFrame(columns=[]).to_csv(comphet_out, sep="\t", index=False)
        print("[INFO] Empty input.", file=sys.stderr)
        return

    header_cols = list(first_chunk.columns)
    if all(col is not None for col in (format_col, proband_col, mother_col, father_col)):
        fmt_c, p_c, m_c, f_c = format_col, proband_col, mother_col, father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(header_cols)

    def make_row_key(row: pd.Series) -> Tuple:
        return (row.get("Chr"), row.get("Start"), row.get("End"), row.get("Ref"), row.get("Alt"))

    def passes_af_filter(series: pd.Series) -> pd.Series:
        af_vals = pd.to_numeric(series, errors="coerce")
        return (af_vals < af_threshold) | (af_vals.isna() if allow_missing_af else False)

    for chunk in [first_chunk] + list(reader):
        if af_col in chunk.columns:
            chunk = chunk[passes_af_filter(chunk[af_col])]
            if chunk.empty:
                continue
        if "Gene.refGene" not in chunk.columns:
            continue

        mask_gene = chunk["Gene.refGene"].apply(
            lambda s: any(g.upper() in genes_of_interest for g in parse_gene_cell(s))
        )
        chunk = chunk[mask_gene]
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

        for i, (fmt, ps, ms, fs) in enumerate(zip(fmt_s, p_s, m_s, f_s)):
            row = chunk.iloc[i]
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                continue

            p_raw = extract_gt_from_sample(ps, gt_idx)
            m_raw = extract_gt_from_sample(ms, gt_idx)
            f_raw = extract_gt_from_sample(fs, gt_idx)

            p_norm = normalize_gt_anyalt(p_raw)
            m_norm = normalize_gt_anyalt(m_raw)
            f_norm = normalize_gt_anyalt(f_raw)
            if p_norm is None:
                continue

            genes = [g for g in parse_gene_cell(row.get("Gene.refGene")) if g.upper() in genes_of_interest]
            if not genes:
                continue

            if p_norm in ("0/1", "1/1"):
                mom_missing = (m_norm is None)
                dad_missing = (f_norm is None)
                both_missing = mom_missing and dad_missing
                if (m_norm == "0/0") and (f_norm == "0/0"):
                    rec = {col: (row[col] if col in row else "") for col in header_cols}
                    rec["dn_category"] = "strict_denovo"
                    denovo_rows.append(rec)
                if both_missing or (mom_missing and f_norm == "0/0") or (dad_missing and m_norm == "0/0"):
                    rec = {col: (row[col] if col in row else "") for col in header_cols}
                    rec["dn_category"] = "dn_missing_parent"
                    denovo_rows.append(rec)
                if (m_norm is not None and f_norm is not None) and is_mendelian_error(p_norm, m_norm, f_norm):
                    rec = {col: (row[col] if col in row else "") for col in header_cols}
                    rec["dn_category"] = "mendelian_error"
                    denovo_rows.append(rec)

            if (p_norm == "0/1") and (m_norm is not None) and (f_norm is not None):
                m_alt = (m_norm in ("0/1", "1/1"))
                f_alt = (f_norm in ("0/1", "1/1"))
                m_is_ref = (m_norm == "0/0")
                f_is_ref = (f_norm == "0/0")
                tag_maternal = m_alt and f_is_ref
                tag_paternal = f_alt and m_is_ref
                if tag_maternal or tag_paternal:
                    rk = make_row_key(row)
                    if rk not in row_store:
                        row_store[rk] = {col: (row[col] if col in row else "") for col in header_cols}
                        row_order[rk] = order_counter
                        order_counter += 1
                    for g in genes:
                        if tag_maternal:
                            gene_maternal.setdefault(g, set()).add(rk)
                        if tag_paternal:
                            gene_paternal.setdefault(g, set()).add(rk)

    if denovo_rows:
        cols = header_cols + (["dn_category"] if "dn_category" not in header_cols else [])
        pd.DataFrame(denovo_rows)[cols].to_csv(denovo_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=header_cols + ["dn_category"]).to_csv(denovo_out, sep="\t", index=False)

    qualifying_genes = sorted(
        g for g in set(gene_maternal) | set(gene_paternal)
        if gene_maternal.get(g) and gene_paternal.get(g)
    )
    rows_to_emit: Set[Tuple] = set()
    for g in qualifying_genes:
        rows_to_emit.update(gene_maternal[g])
        rows_to_emit.update(gene_paternal[g])

    if rows_to_emit:
        out_rows = [row_store[rk] for rk in sorted(rows_to_emit, key=lambda k: row_order.get(k, 10**18))]
        pd.DataFrame(out_rows, columns=header_cols).to_csv(comphet_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=header_cols).to_csv(comphet_out, sep="\t", index=False)

    print("Scan complete.")
    print(f"  De novo/M.E. candidates: {len(denovo_rows)} -> {denovo_out}")
    print(f"  Genes with maternal-only tags: {sum(1 for g,s in gene_maternal.items() if s and not gene_paternal.get(g))}")
    print(f"  Genes with paternal-only tags: {sum(1 for g,s in gene_paternal.items() if s and not gene_maternal.get(g))}")
    print(f"  Genes meeting comp-het (both sides): {len(qualifying_genes)}")
    print(f"  Comp-het variants emitted: {len(rows_to_emit)} -> {comphet_out}")


def main():
    ap = argparse.ArgumentParser(
        description="Find candidate de novo (incl. ME & missing-parent) and compound-heterozygous variants in genes of interest from ANNOVAR trio TSV."
    )
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output-prefix", required=True)
    ap.add_argument("--genes")
    ap.add_argument("--format-col", default=None)
    ap.add_argument("--proband-col", default=None)
    ap.add_argument("--mother-col", default=None)
    ap.add_argument("--father-col", default=None)
    ap.add_argument("--af-col", default="gnomad41_genome_AF")
    ap.add_argument("--af-threshold", type=float, default=0.001)
    ap.add_argument("--require-af", action="store_true")
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()

    genes = load_gene_set(args.genes)
    process_trio(
        input_file=args.input,
        output_prefix=args.output_prefix,
        genes_of_interest=genes,
        format_col=args.format_col,
        proband_col=args.proband_col,
        mother_col=args.mother_col,
        father_col=args.father_col,
        af_col=args.af_col,
        af_threshold=args.af_threshold,
        allow_missing_af=(not args.require_af),
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()

