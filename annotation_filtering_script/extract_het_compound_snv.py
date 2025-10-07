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


def exonic_consequence_ok(exonic_func: str) -> bool:
    if not isinstance(exonic_func, str):
        return False
    s = exonic_func.strip().lower()
    return s not in ("", ".", "synonymous snv", "unknown")


def is_exonic_func_refgene(func_refgene: str) -> bool:
    if not isinstance(func_refgene, str):
        return False
    return "exonic" in func_refgene.strip().lower()


def is_sex_chrom(chrom: str) -> bool:
    if not isinstance(chrom, str):
        return False
    c = chrom.strip()
    c = c[3:] if c.lower().startswith("chr") else c
    return c.upper() in ("X", "Y")


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


def process_comphet_lenient(
    input_file: str,
    output_prefix: str,
    format_col: Optional[str] = None,
    proband_col: Optional[str] = None,
    mother_col: Optional[str] = None,
    father_col: Optional[str] = None,
    af_col: str = "gnomad41_genome_AF",
    af_threshold: float = 0.001,
    chunksize: int = 100000,
):
    out_file = f"{output_prefix}_comphet_lenient.tsv"
    header_cols: List[str] = []
    gt_index_cache: Dict[str, int] = {}
    gene_maternal: Dict[str, Set[Tuple]] = {}
    gene_paternal: Dict[str, Set[Tuple]] = {}
    row_store: Dict[Tuple, Dict[str, str]] = {}
    row_order: Dict[Tuple, int] = {}
    order_counter = 0
    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)
    try:
        first_chunk = next(reader)
    except StopIteration:
        print("[INFO] Empty input.", file=sys.stderr)
        return
    header_cols = list(first_chunk.columns)
    if all(col is not None for col in (format_col, proband_col, mother_col, father_col)):
        fmt_c, p_c, m_c, f_c = format_col, proband_col, mother_col, father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(header_cols)

    def make_row_key(row: pd.Series) -> Tuple:
        return (row.get("Chr"), row.get("Start"), row.get("End"), row.get("Ref"), row.get("Alt"))

    for chunk in [first_chunk] + list(reader):
        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        chunk = chunk[af_vals.isna() | (af_vals < af_threshold)]
        if chunk.empty:
            continue
        if "Func.refGene" in chunk.columns:
            chunk = chunk[chunk["Func.refGene"].apply(is_exonic_func_refgene)]
            if chunk.empty:
                continue
        if "ExonicFunc.refGene" in chunk.columns:
            chunk = chunk[chunk["ExonicFunc.refGene"].apply(exonic_consequence_ok)]
            if chunk.empty:
                continue
        for col in (fmt_c, p_c, m_c, f_c, "Chr", "Gene.refGene"):
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
            if is_sex_chrom(row["Chr"]):
                continue
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                continue
            p_norm = normalize_gt_anyalt(extract_gt_from_sample(ps, gt_idx))
            m_norm = normalize_gt_anyalt(extract_gt_from_sample(ms, gt_idx))
            f_norm = normalize_gt_anyalt(extract_gt_from_sample(fs, gt_idx))
            if p_norm != "0/1":
                continue
            m_alt = (m_norm in ("0/1", "1/1"))
            f_alt = (f_norm in ("0/1", "1/1"))
            m_noncarrier = (m_norm in ("0/0", None))
            f_noncarrier = (f_norm in ("0/0", None))
            tag_maternal = m_alt and f_noncarrier
            tag_paternal = f_alt and m_noncarrier
            if not (tag_maternal or tag_paternal):
                continue
            genes_cell = row.get("Gene.refGene")
            if not isinstance(genes_cell, str) or genes_cell.strip() in ("", "."):
                continue
            genes = [g.strip() for g in genes_cell.split(";") if g.strip() not in ("", ".")]
            if not genes:
                continue
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

    qualifying_genes = {g for g in gene_maternal if gene_maternal[g] and gene_paternal.get(g)}
    rows_to_emit: Set[Tuple] = set()
    for g in qualifying_genes:
        rows_to_emit.update(gene_maternal[g])
        rows_to_emit.update(gene_paternal[g])

    if rows_to_emit:
        out_rows = [row_store[rk] for rk in sorted(rows_to_emit, key=lambda k: row_order.get(k, 10**18))]
        pd.DataFrame(out_rows, columns=header_cols).to_csv(out_file, sep="\t", index=False)
    else:
        pd.DataFrame(columns=header_cols).to_csv(out_file, sep="\t", index=False)

    print("Compound-het (lenient) scan complete.")
    print(f"  Genes with maternal-only candidates: {sum(1 for g,s in gene_maternal.items() if s and not gene_paternal.get(g))}")
    print(f"  Genes with paternal-only candidates: {sum(1 for g,s in gene_paternal.items() if s and not gene_maternal.get(g))}")
    print(f"  Genes meeting lenient comp-het (both sides): {len(qualifying_genes)}")
    print(f"  Variants emitted: {len(rows_to_emit)} -> {out_file}")


def main():
    ap = argparse.ArgumentParser(description="Identify inherited compound-heterozygous candidates (lenient) from ANNOVAR trio TSV.")
    ap.add_argument("-i", "--input",
    ap.add_argument("-o", "--output-prefix", required=True)
    ap.add_argument("--format-col", default=None)
    ap.add_argument("--proband-col", default=None)
    ap.add_argument("--mother-col", default=None)
    ap.add_argument("--father-col", default=None)
    ap.add_argument("--af-col", default="gnomad41_genome_AF")
    ap.add_argument("--af-threshold", type=float, default=0.001)
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()

    process_comphet_lenient(
        input_file=args.input,
        output_prefix=args.output_prefix,
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
